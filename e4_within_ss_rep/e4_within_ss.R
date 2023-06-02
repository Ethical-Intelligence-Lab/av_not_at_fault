## ================================================================================================================
##                                 Harvard Business School, Ethical Intelligence Lab
## ================================================================================================================
##                                DATA ANALYSIS | AV SCENARIOS | WITHIN SUBJECTS               
## ================================================================================================================

## clear workspace
rm(list = ls()) 

options(download.file.method="libcurl")

## install packages
library(ggpubr)
library(rstatix)
if (!require(pacman)) {install.packages("pacman")}
pacman::p_load('ggplot2',         # plotting
               'ggsignif',        # plotting significance bars
               'lme4',            # functions for fitting linear regression models
               'ggforce',         # make ggplot even fancier
               'ggpubr',          # arrange plots in a grid, if needed
               'ltm',             # probably not using..
               'tidyr',           # tools for cleaning messy data
               'stringr',         # perform string substitutions easily
               'assertthat',      # allows me to check whether a variable is a string, with is.string
               'lsmeans',         # contrast analysis for regression models
               'stats',           # use function to adjust for multiple comparisons
               'filesstrings',    # create and move files
               'simr',            # power analysis for mixed models
               'compute.es',      # effect size package
               'effsize',         # another effect size package
               'pwr',             # package for power calculation
               'nlme',            # get p values for mixed effect model
               'DescTools',       # get Cramer's V
               'Hmisc'
)

library("lmerTest")

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
# if importing from Qualtrics: (i) export data as numeric values, and (ii) delete rows 2 and 3 of the .csv file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
d <- read.csv('e4_within_ss_50.csv') 

## explore dataframe: 
dim(d) # will provide dimensions of the dataframe by row [1] and column [2]
colnames(d) # will provide all column names
summary(d)
dim(d)[1]

## rename condition variables:
# names(d)[names(d) == 'order'] <- 'agent_order'

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d)

## ================================================================================================================
##                                              PERFORM EXCLUSIONS                
## ================================================================================================================

## get number of participants BEFORE exclusions: 
n_original <- dim(d)[1] # extracting number of rows only, not columns
n_original 

## perform comprehension exclusions: 
# this will remove responses from the dataframe that failed comprehension checks (i.e., "2")
d <- subset(d, ( d$comp_1_AV == 1 & d$comp_2_AV == 1 & d$comp_1_HDV == 2 & d$comp_2_HDV == 1))
dim(d) # number of participants should decrease after comprehension exclusions

## get number of participants AFTER exclusions: 
n_final <- dim(d)[1] # extracting number of rows only, not columns
n_final 
percent_excluded <- (n_original - n_final)/n_original 
percent_excluded

## duplicate AV condition vB liability column to match with HDV condition driver liability col
d$vB_mnfctr_liable_AV_2 <- d$vB_mnfctr_liable_AV_1
d <- d %>% relocate(vB_mnfctr_liable_AV_2, .after=vB_mnfctr_liable_AV_1)

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_subset <- array(dim=c(dim(d)[1], 13))
colnames(d_subset) <- c('order', 
                        'vA_liable_AV', 'vB_m_v_d_liable_AV', 'vB_m_v_m_liable_AV',
                        'vA_cntrfctl_AV', 'vB_cntrfctl_AV', 'avoid_AV',
                        'vA_liable_HDV', 'vB_m_v_d_liable_HDV', 'vB_m_v_m_liable_HDV',
                        'vA_cntrfctl_HDV', 'vB_cntrfctl_HDV', 'avoid_HDV')
d_subset <- as.data.frame(d_subset, stringsAsFactors=FALSE) 

## extract good data from the middle part of raw data:
for(i in 1:dim(d)[1]) {
    curr_AV <- d[i,c(29:34)][!is.na(d[i,c(29:34)])] # for a given row, get only the non-NA values
    curr_HDV <- d[i,c(45:51)][!is.na(d[i,c(45:51)])]
    d_subset[i,2:7] <- as.numeric(curr_AV[curr_AV!= ""])# and only the non-empty values
    d_subset[i,8:13] <- as.numeric(curr_HDV[curr_HDV!= ""])
    d_subset[i,1] <- d[i,'order']
}

d_merged <- array(dim=c(dim(d)[1]*2, 8))
colnames(d_merged) <- c('order', 'cond',
                        'vA_liable', 'vB_m_v_d_liable', 'vB_m_v_m_liable',
                        'vA_cntrfctl', 'vB_cntrfctl', 'avoid')
d_merged <- as.data.frame(d_merged, stringsAsFactors=FALSE) 
d_merged$order <- c(d_subset$order, d_subset$order)
d_merged$cond <- c(rep(c(1),each=dim(d)[1]), rep(c(2),each=dim(d)[1])) #1= AI, 2 =Human

d_merged$vA_liable <- c(d_subset$vA_liable_AV, d_subset$vA_liable_HDV)
d_merged$vB_m_v_d_liable <- c(d_subset$vB_m_v_d_liable_AV, d_subset$vB_m_v_d_liable_HDV)
d_merged$vB_m_v_m_liable <- c(d_subset$vB_m_v_m_liable_AV, d_subset$vB_m_v_m_liable_HDV)
d_merged$vA_cntrfctl <- c(d_subset$vA_cntrfctl_AV, d_subset$vA_cntrfctl_HDV)
d_merged$vB_cntrfctl <- c(d_subset$vB_cntrfctl_AV, d_subset$vB_cntrfctl_HDV)
d_merged$avoid <- c(d_subset$avoid_AV, d_subset$avoid_HDV)

# creates a new variable for the string version of condition
d_merged$cond_name <- ifelse(d_merged$cond==1, "av", "human")

## ================================================================================================================
##                                            PARTICIPANT CHARACTERISTICS                 
## ================================================================================================================

# ## % prior experience apps
# table(d$ai_familiarity)[1]/sum(table(d$ai_familiarity))
# d$ai_familiarity
# #more like 11
# 11/dim(d)[1]

## age
mean(as.numeric(d$age), trim = 0, na.rm = TRUE) ## mean age 

## gender
table(d$gender)[2]/sum(table(d$gender)) ## percentage of females

## ================================================================================================================
##                              DATA ANALYSIS - SUMMARY, T-TESTS, AND LINEAR REGRESSION               
## ================================================================================================================

## (1) LIABLE VEHICLE A DRIVER
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vA_liable, type = "mean_sd")
mean(as.numeric(d_subset$vA_liable_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vA_liable_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vA_liable_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vA_liable_HDV), na.rm = TRUE)

vA_liable_T <- t.test(vA_liable ~ cond, data = d_merged, paired = TRUE) 
vA_liable_T

vA_liable_mod <- lmer(vA_liable ~ cond + (1 | order), data=d_merged)
summary(vA_liable_mod)

## (2) LIABLE VEHICLE B MANUFACTURER VS LIABLE HDV DRIVER
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vB_m_v_d_liable, type = "mean_sd")
mean(as.numeric(d_subset$vB_m_v_d_liable_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_m_v_d_liable_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vB_m_v_d_liable_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_m_v_d_liable_HDV), na.rm = TRUE)

vB_m_v_d_liable_T <- t.test(vB_m_v_d_liable ~ cond, data = d_merged, paired = TRUE) 
vB_m_v_d_liable_T

vB_m_v_d_liable_mod <- lmer(vB_m_v_d_liable ~ cond + (1 | order), data=d_merged)
summary(vB_m_v_d_liable_mod)

## (3) LIABLE VEHICLE B MANUFACTURER VS LIABLE HDV MANUFACTURER
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vB_m_v_m_liable, type = "mean_sd")
mean(as.numeric(d_subset$vB_m_v_m_liable_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_m_v_m_liable_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vB_m_v_m_liable_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_m_v_m_liable_HDV), na.rm = TRUE)

vB_m_v_m_liable_T <- t.test(vB_m_v_m_liable ~ cond, data = d_merged, paired = TRUE) 
vB_m_v_m_liable_T

vB_m_v_m_liable_mod <- lmer(vB_m_v_m_liable ~ cond + (1 | order), data=d_merged)
summary(vB_m_v_m_liable_mod)

## (4) CONSIDER VEHICLE A COUNTERFACTUAL
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vA_cntrfctl, type = "mean_sd")
mean(as.numeric(d_subset$vA_cntrfctl_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vA_cntrfctl_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vA_cntrfctl_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vA_cntrfctl_HDV), na.rm = TRUE)

vA_cntrfctl_T <- t.test(vA_cntrfctl ~ cond, data = d_merged, paired = TRUE) 
vA_cntrfctl_T

vA_cntrfctl_mod <- lmer(vA_cntrfctl ~ cond + (1 | order), data=d_merged)
summary(vA_cntrfctl_mod)

## (5) CONSIDER VEHICLE B COUNTERFACTUAL
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(vB_cntrfctl, type = "mean_sd")
mean(as.numeric(d_subset$vB_cntrfctl_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_cntrfctl_AV), na.rm = TRUE)
mean(as.numeric(d_subset$vB_cntrfctl_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$vB_cntrfctl_HDV), na.rm = TRUE)

vB_cntrfctl_T <- t.test(vB_cntrfctl ~ cond, data = d_merged, paired = TRUE) 
vB_cntrfctl_T

vB_cntrfctl_mod <- lmer(vB_cntrfctl ~ cond + (1 | order), data=d_merged)
summary(vB_cntrfctl_mod)

## (6) VEHICLE B CAN AVOID
## get summary statistics
d_merged %>% group_by(cond_name) %>% get_summary_stats(avoid, type = "mean_sd")
mean(as.numeric(d_subset$avoid_AV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$avoid_AV), na.rm = TRUE)
mean(as.numeric(d_subset$avoid_HDV), trim = 0, na.rm = TRUE)
sd(as.numeric(d_subset$avoid_HDV), na.rm = TRUE)

avoid_T <- t.test(avoid ~ cond, data = d_merged, paired = TRUE) 
avoid_T

avoid_mod <- lmer(avoid ~ cond + (1 | order), data=d_merged)
summary(avoid_mod)

## ================================================================================================================
##                                              MEDIATION ANALYSIS                
## ================================================================================================================

source("../process.R")

# SERIAL MEDIATION
process(data = d_merged, y = "vB_m_v_m_liable", x = "cond", 
        m =c("vB_cntrfctl", "avoid"), model = 6, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
process(data = d_merged, y = "vB_m_v_d_liable", x = "cond", 
        m =c("vB_cntrfctl", "avoid"), model = 6, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)


## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
t_names <- c("AV", "HDV")
title_size <- 20

## function for getting the correct sig annotation
get_annotation <- function(p_val) {
  if (p_val < 0.001) {
    return (list('***', 5.5))
  } else if (p_val < 0.01) {
    return (list('**', 5.5))
  } else if (p_val < 0.05) {
    return (list('*', 5.5))
  } else if (p_val < 0.1) {
    return (list('^', 5.5))
  } else {
    return (list('NS', 3))
  }
}

## (1) VA driver liable
annotations <- get_annotation(vA_liable_T$p.value)
p1_1 <- ggplot(d_merged,aes(x=factor(cond_name),y=vA_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_1 <- p1_1 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. A Driver Liability") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_1

## (2) VB manufacturer/driver liability
annotations <- get_annotation(vB_m_v_d_liable_T$p.value)
p1_2 <- ggplot(d_merged,aes(x=factor(cond_name),y=vB_m_v_d_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_2 <- p1_2 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Manufacturer\nor Driver Liability") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  #theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_2

## (3) VB manufacturer liability
annotations <- get_annotation(vB_m_v_m_liable_T$p.value)
p1_3 <- ggplot(d_merged,aes(x=factor(cond_name),y=vB_m_v_m_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_3 <- p1_3 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Manufacturer Liability") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_3

## (4) VA counterfactual
annotations <- get_annotation(vA_cntrfctl_T$p.value)
p1_4 <- ggplot(d_merged,aes(x=factor(cond_name),y=vA_cntrfctl)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_4 <- p1_4 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Consider Veh. A Counterfactual") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_4

## (5) VB Counterfactual
annotations <- get_annotation(vB_cntrfctl_T$p.value)
p1_5 <- ggplot(d_merged,aes(x=factor(cond_name),y=vB_cntrfctl)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_5 <- p1_5 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Consider Veh. B Counterfactual") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_5

## (6) Capability to Avoid
annotations <- get_annotation(avoid_T$p.value)
p1_6 <- ggplot(d_merged,aes(x=factor(cond_name),y=avoid)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("av", "human")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_6 <- p1_6 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Capability to Avoid") +
  xlab ("") + ylab ("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(plot.title = element_text(size=12, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_6

## PLOT SERIES 1
dev.new(width=10,height=8,noRStudioGD = TRUE)
figure1 <- ggarrange(p1_1, p1_2, p1_3, p1_4, p1_5, p1_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure1 <- annotate_figure(figure1,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 

plot(figure1)

write.csv(d_merged, 'd_spss.csv')

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================