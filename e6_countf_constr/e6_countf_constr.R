## clear workspace
rm(list = ls()) 

source("../process.R")

options(download.file.method="libcurl")

## install packages
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
               'dplyr',           # package to move columns around
               'Hmisc'
)

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
# if importing from Qualtrics: (i) export data as numeric values, and (ii) delete rows 2 and 3 of the .csv file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
d <- read.csv('e6_countf_constr_300.csv')

## explore dataframe: 
dim(d) # will provide dimensions of the dataframe by row [1] and column [2]
colnames(d) # will provide all column names
summary(d)

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d) # number of participants should decrease after attention exclusions

## split up dataframes between AV and HDV conditions
## this is necessary before comprehension exclusions
d_AV_cnstr <- subset(d, (d$FL_4_DO == "FL_39"))
d_HDV_cnstr <- subset(d, (d$FL_4_DO == "FL_40"))
d_AV_uncnstr <- subset(d, (d$FL_4_DO == "FL_54"))
d_HDV_uncnstr <- subset(d, (d$FL_4_DO == "FL_58"))

## get number of participants BEFORE exclusions: 
n_original <- dim(d)[1] # extracting number of rows only, not columns
n_original_AV_cnstr <- dim(d_AV_cnstr)[1]
n_original_HDV_cnstr <- dim(d_HDV_cnstr)[1]
n_original_AV_uncnstr <- dim(d_AV_uncnstr)[1]
n_original_HDV_uncnstr <- dim(d_HDV_uncnstr)[1]

## perform comprehension exclusions separately for AV and HDV: 
# this will remove responses from the dataframe that failed comprehension checks (i.e., "2")
d_AV_cnstr <- subset(d_AV_cnstr, (d_AV_cnstr$comp1 == 1 & d_AV_cnstr$comp_accident == 1))
d_HDV_cnstr <- subset(d_HDV_cnstr, (d_HDV_cnstr$comp1 == 2 & d_HDV_cnstr$comp_accident == 1))
d_AV_uncnstr <- subset(d_AV_uncnstr, (d_AV_uncnstr$comp1 == 1 & d_AV_uncnstr$comp_accident == 1))
d_HDV_uncnstr <- subset(d_HDV_uncnstr, (d_HDV_uncnstr$comp1 == 2 & d_HDV_uncnstr$comp_accident == 1))
dim(d_AV_cnstr) # number of participants should decrease after comprehension exclusions
dim(d_HDV_cnstr)
dim(d_AV_uncnstr)
dim(d_HDV_uncnstr)

## get number of participants AFTER exclusions: 
n_final_AV_cnstr <- dim(d_AV_cnstr)[1] # extracting number of rows only, not columns
n_final_HDV_cnstr <- dim(d_HDV_cnstr)[1]
n_final_AV_uncnstr <- dim(d_AV_uncnstr)[1]
n_final_HDV_uncnstr <- dim(d_HDV_uncnstr)[1]
n_final <- n_final_AV_cnstr + n_final_HDV_cnstr + n_final_AV_uncnstr + n_final_HDV_uncnstr
percent_excluded <- (n_original - n_final)/n_original
percent_excluded_AV_cnstr <- (n_original_AV_cnstr - n_final_AV_cnstr)/n_original_AV_cnstr 
percent_excluded_HDV_cnstr <- (n_original_HDV_cnstr - n_final_HDV_cnstr)/n_original_HDV_cnstr
percent_excluded_AV_uncnstr <- (n_original_AV_uncnstr - n_final_AV_uncnstr)/n_original_AV_uncnstr 
percent_excluded_HDV_uncnstr <- (n_original_HDV_uncnstr - n_final_HDV_uncnstr)/n_original_HDV_uncnstr

## remove unused columns (other condition and click info) according to condition
d_AV_cnstr <- d_AV_cnstr[-c(21:28,34:74)] # first=click info, second=all other columns
d_HDV_cnstr <- d_HDV_cnstr[-c(21:33,34:41,48:74)] # first=prior columns, second=click info, third=later columns
d_AV_uncnstr <- d_AV_uncnstr[-c(21:47,48:55,61:74)] # first=prior columns, second=click info, third=later columns
d_HDV_uncnstr <- d_HDV_uncnstr[-c(21:60,61:68)] # first=prior columns, second=click info


## duplicate AV condition vB liability column to match with HDV condition driver liability col
d_AV_cnstr$vB_frm_liab_AV_cnstr_2 <- d_AV_cnstr$vB_frm_liab_AV_cnstr_1
d_AV_cnstr <- d_AV_cnstr %>% relocate(vB_frm_liab_AV_cnstr_2, .after=vB_frm_liab_AV_cnstr_1)
d_AV_uncnstr$vB_frm_liab_AV_2 <- d_AV_uncnstr$vB_frm_liab_AV_1
d_AV_uncnstr <- d_AV_uncnstr %>% relocate(vB_frm_liab_AV_2, .after=vB_frm_liab_AV_1)

## get mean age and gender:
mean_age = mean(as.numeric(d$age), na.rm = TRUE) # removing NAs from the dataframe before computing mean 
gender = table(d$gender)[1]/sum(table(d$gender)) # percent male

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_subset <- array(dim=c(n_final, 10))
colnames(d_subset) <- c('agent_name', 'scen_name', 'vA_liable', 'vB_m_v_d_liable', 'vB_m_v_m_liable', 'vA_cntrfctl', 'vB_cntrfctl', 
                        'avoid', 'comp1', 'comp2')
d_subset <- as.data.frame(d_subset, stringsAsFactors=FALSE) 

## extract good data from the middle part of raw data in AV constrained
for(i in 1:n_final_AV_cnstr) {
  curr <- d_AV_cnstr[i,21:28][!is.na(d_AV_cnstr[i,21:28])] # for a given row, get only the non-NA values
  d_subset[i,3:10] <- as.numeric(curr[curr!= ""]) # and only the non-empty values
  #d_subset[i,1] <- d_AV[i,43][!is.na(d_AV[i,43])]
  d_subset[i,1] <- "av"
  d_subset[i,2] <- "cnstr"
}

## extract good data from the middle part of raw data in HDV constrained
for(i in 1:n_final_HDV_cnstr) {
  j = i+n_final_AV_cnstr
  curr <- d_HDV_cnstr[i,21:28][!is.na(d_HDV_cnstr[i,21:28])] # for a given row, get only the non-NA values
  d_subset[j,3:10] <- as.numeric(curr) # and only the non-empty values
  d_subset[j,1] <- "human"
  d_subset[j,2] <- "cnstr"
}

## extract good data from the middle part of raw data in AV unconstrained
for(i in 1:n_final_AV_uncnstr) {
  j = i+n_final_AV_cnstr+n_final_HDV_cnstr
  curr <- d_AV_uncnstr[i,21:28][!is.na(d_AV_uncnstr[i,21:28])] # for a given row, get only the non-NA values
  d_subset[j,3:10] <- as.numeric(curr) # and only the non-empty values
  d_subset[j,1] <- "av"
  d_subset[j,2] <- "uncnstr"
}

## extract good data from the middle part of raw data in HDV unconstrained
for(i in 1:n_final_HDV_uncnstr) {
  j = i+n_final_AV_cnstr+n_final_HDV_cnstr+n_final_AV_uncnstr
  curr <- d_HDV_uncnstr[i,21:28][!is.na(d_HDV_uncnstr[i,21:28])] # for a given row, get only the non-NA values
  d_subset[j,3:10] <- as.numeric(curr) # and only the non-empty values
  d_subset[j,1] <- "human"
  d_subset[j,2] <- "uncnstr"
}

## just to keep the df names straight for next section
d_merged <- d_subset

# cond_n where av=1, human=2
d_merged$agent_n <- ifelse(d_merged$agent_name=="av", 1, 2)
d_merged$scen_n <- ifelse(d_merged$scen_name=="cnstr", 1, 2)

## ================================================================================================================
##                                             DATA ANALYSIS - T-TESTS               
## ================================================================================================================

## MANUFACTURER VS DRIVER --------
## get summary statistics
d_merged %>%
  group_by(scen_n) %>%
  get_summary_stats(vB_m_v_d_liable, type = "mean_sd")

## anova
m_v_d_mod <- aov(vB_m_v_d_liable ~ as.factor(agent_n) * as.factor(scen_n), data = d_merged)
summary(m_v_d_mod)

## MANUFACTURER VS MANUFACTURER --------
## get summary statistics
d_merged %>%
  group_by(scen_n) %>%
  get_summary_stats(vB_m_v_m_liable, type = "mean_sd")

## anova
m_v_m_mod <- aov(vB_m_v_m_liable ~ as.factor(agent_n) * as.factor(scen_n), data = d_merged)
summary(m_v_m_mod)

table(d_merged$con) #give us table of number of people in each condition - want to have equal number of people in each condition

## (1) LIABLE VEHICLE A DRIVER
vA_liable_T_agent <- t.test(vA_liable ~ agent_name, data = d_merged, paired = FALSE) 
vA_liable_T_agent$parameter
vA_liable_T_agent$statistic
vA_liable_T_agent$p.value

vA_liable_T_scen <- t.test(vA_liable ~ scen_name, data = d_merged, paired = FALSE) 
vA_liable_T_scen$parameter
vA_liable_T_scen$statistic
vA_liable_T_scen$p.value

vA_liable_T_cnstr <- t.test(d_merged$vA_liable[d_merged$scen_name=="cnstr" & d_merged$agent_name == "av"], 
                            d_merged$vA_liable[d_merged$scen_name=="cnstr" & d_merged$agent_name == "human"], paired=FALSE)
vA_liable_T_uncnstr <- t.test(d_merged$vA_liable[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "av"], 
                              d_merged$vA_liable[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "human"], paired=FALSE)

## (2) LIABLE VEHICLE B MANUFACTURER VS LIABLE HDV DRIVER
vB_m_v_d_liable_T_agent <- t.test(vB_m_v_d_liable ~ agent_name, data = d_merged, paired = FALSE) 
vB_m_v_d_liable_T_agent$parameter
vB_m_v_d_liable_T_agent$statistic
vB_m_v_d_liable_T_agent$p.value

vB_m_v_d_liable_T_scen <- t.test(vB_m_v_d_liable ~ scen_name, data = d_merged, paired = FALSE) 
vB_m_v_d_liable_T_scen$parameter
vB_m_v_d_liable_T_scen$statistic
vB_m_v_d_liable_T_scen$p.value

vB_m_v_d_liable_T_cnstr <- t.test(d_merged$vB_m_v_d_liable[d_merged$scen_name=="cnstr" & d_merged$agent_name == "av"], 
                                  d_merged$vB_m_v_d_liable[d_merged$scen_name=="cnstr" & d_merged$agent_name == "human"], paired=FALSE)
vB_m_v_d_liable_T_uncnstr <- t.test(d_merged$vB_m_v_d_liable[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "av"], 
                                    d_merged$vB_m_v_d_liable[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "human"], paired=FALSE)

## (3) LIABLE VEHICLE B MANUFACTURER VS LIABLE HDV MANUFACTURER
vB_m_v_m_liable_T_agent <- t.test(vB_m_v_m_liable ~ agent_name, data = d_merged, paired = FALSE) 
vB_m_v_m_liable_T_agent$parameter
vB_m_v_m_liable_T_agent$statistic
vB_m_v_m_liable_T_agent$p.value

vB_m_v_m_liable_T_scen <- t.test(vB_m_v_m_liable ~ scen_name, data = d_merged, paired = FALSE) 
vB_m_v_m_liable_T_scen$parameter
vB_m_v_m_liable_T_scen$statistic
vB_m_v_m_liable_T_scen$p.value

vB_m_v_m_liable_T_cnstr <- t.test(d_merged$vB_m_v_m_liable[d_merged$scen_name=="cnstr" & d_merged$agent_name == "av"], 
                                  d_merged$vB_m_v_m_liable[d_merged$scen_name=="cnstr" & d_merged$agent_name == "human"], paired=FALSE)
vB_m_v_m_liable_T_uncnstr <- t.test(d_merged$vB_m_v_m_liable[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "av"], 
                                    d_merged$vB_m_v_m_liable[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "human"], paired=FALSE)

## (4) CONSIDER VEHICLE A COUNTERFACTUAL
vA_cntrfctl_T_agent <- t.test(vA_cntrfctl ~ agent_name, data = d_merged, paired = FALSE) 
vA_cntrfctl_T_agent$parameter
vA_cntrfctl_T_agent$statistic
vA_cntrfctl_T_agent$p.value

vA_cntrfctl_T_scen <- t.test(vA_cntrfctl ~ scen_name, data = d_merged, paired = FALSE) 
vA_cntrfctl_T_scen$parameter
vA_cntrfctl_T_scen$statistic
vA_cntrfctl_T_scen$p.value

vA_cntrfctl_T_cnstr <- t.test(d_merged$vA_cntrfctl[d_merged$scen_name=="cnstr" & d_merged$agent_name == "av"], 
                              d_merged$vA_cntrfctl[d_merged$scen_name=="cnstr" & d_merged$agent_name == "human"], paired=FALSE)
vA_cntrfctl_T_uncnstr <- t.test(d_merged$vA_cntrfctl[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "av"], 
                                d_merged$vA_cntrfctl[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "human"], paired=FALSE)

## (5) CONSIDER VEHICLE B COUNTERFACTUAL
vB_cntrfctl_T_agent <- t.test(vB_cntrfctl ~ agent_name, data = d_merged, paired = FALSE) 
vB_cntrfctl_T_agent$parameter
vB_cntrfctl_T_agent$statistic
vB_cntrfctl_T_agent$p.value

vB_cntrfctl_T_scen <- t.test(vB_cntrfctl ~ scen_name, data = d_merged, paired = FALSE) 
vB_cntrfctl_T_scen$parameter
vB_cntrfctl_T_scen$statistic
vB_cntrfctl_T_scen$p.value

vB_cntrfctl_T_cnstr <- t.test(d_merged$vB_cntrfctl[d_merged$scen_name=="cnstr" & d_merged$agent_name == "av"], 
                              d_merged$vB_cntrfctl[d_merged$scen_name=="cnstr" & d_merged$agent_name == "human"], paired=FALSE)
vB_cntrfctl_T_uncnstr <- t.test(d_merged$vB_cntrfctl[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "av"], 
                                d_merged$vB_cntrfctl[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "human"], paired=FALSE)

## (6) VEHICLE B CAN AVOID
avoid_T_agent <- t.test(avoid ~ agent_name, data = d_merged, paired = FALSE) 
avoid_T_agent$parameter
avoid_T_agent$statistic
avoid_T_agent$p.value

avoid_T_scen <- t.test(avoid ~ scen_name, data = d_merged, paired = FALSE) 
avoid_T_scen$parameter
avoid_T_scen$statistic
avoid_T_scen$p.value

avoid_T_cnstr <- t.test(d_merged$avoid[d_merged$scen_name=="cnstr" & d_merged$agent_name == "av"], 
                        d_merged$avoid[d_merged$scen_name=="cnstr" & d_merged$agent_name == "human"], paired=FALSE)
avoid_T_uncnstr <- t.test(d_merged$avoid[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "av"], 
                        d_merged$avoid[d_merged$scen_name=="uncnstr" & d_merged$agent_name == "human"], paired=FALSE)

cor(d_merged[,3:8])

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

# MODERATED MEDIATION
#  the effect of agent on judgments of liability
process(data = d_merged, y = "vB_m_v_m_liable", x = "agent_n",
        m =c("vB_cntrfctl"), w = "scen_n", model = 7, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
process(data = d_merged, y = "vB_m_v_d_liable", x = "agent_n",
        m =c("vB_cntrfctl"), w = "scen_n", model = 7, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

# MODERATED SERIAL MEDIATION
# 87 = B path, 83 = A path, 91 = center path
process(data = d_merged, y = "vB_m_v_m_liable", x = "agent_n", 
        m =c("vB_cntrfctl", "avoid"), w = "scen_n", model = 83, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
process(data = d_merged, y = "vB_m_v_d_liable", x = "agent_n", 
        m =c("vB_cntrfctl", "avoid"), w = "scen_n", model = 83, effsize =1, total =1, stand =1, 
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
## FL39 --> AV condition; FL40 --> HDV condition
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
annotations <- get_annotation(vA_liable_T_agent$p.value)
p1_1 <- ggplot(d_merged,aes(x=factor(agent_name),y=vA_liable)) +  
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
annotations <- get_annotation(vB_m_v_d_liable_T_agent$p.value)
p1_2 <- ggplot(d_merged,aes(x=factor(agent_name),y=vB_m_v_d_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncsntr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

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
annotations <- get_annotation(vB_m_v_m_liable_T_agent$p.value)
p1_3 <- ggplot(d_merged,aes(x=factor(agent_name),y=vB_m_v_m_liable)) +  
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
annotations <- get_annotation(vA_cntrfctl_T_agent$p.value)
p1_4 <- ggplot(d_merged,aes(x=factor(agent_name),y=vA_cntrfctl)) +  
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
annotations <- get_annotation(vB_cntrfctl_T_agent$p.value)
p1_5 <- ggplot(d_merged,aes(x=factor(agent_name),y=vB_cntrfctl)) +  
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
annotations <- get_annotation(avoid_T_agent$p.value)
p1_6 <- ggplot(d_merged,aes(x=factor(agent_name),y=avoid)) +  
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

## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
t_names <- c("Constrained", "Unconstrained")

# (1) VA driver liable
cnstr_anno <- get_annotation(vA_liable_T_cnstr$p.value)
uncnstr_anno <- get_annotation(vA_liable_T_uncnstr$p.value)
dev.new(width=13,height=6,noRStudioGD = TRUE)
p3_1 <- ggplot(d_merged,aes(x=factor(scen_name),y=vA_liable, fill=agent_name)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(cnstr_anno[1]),unlist(uncnstr_anno[1])), textsize=7.5)
p3_1 <- p3_1 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. A Driver Liability") +
  xlab ("Scenario Type") + ylab ("Measure") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=18, hjust=0.5)) +
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14), legend.position="top")+
  labs(fill='')+
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3_1

## (2) VB manufacturer/driver liability
cnstr_anno <- get_annotation(vB_m_v_d_liable_T_cnstr$p.value)
uncnstr_anno <- get_annotation(vB_m_v_d_liable_T_uncnstr$p.value)
dev.new(width=13,height=6,noRStudioGD = TRUE)
p3_2 <- ggplot(d_merged,aes(x=factor(scen_name),y=vB_m_v_d_liable, fill=agent_name)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(cnstr_anno[1]),unlist(uncnstr_anno[1])), textsize=7.5)
p3_2 <- p3_2 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Manufacturer\nor Driver Liability") +
  xlab ("Scenario Type") + ylab ("Measure") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=18, hjust=0.5)) +
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14), legend.position="top")+
  labs(fill='')+
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3_2

## (3) VB manufacturer liability
cnstr_anno <- get_annotation(vB_m_v_m_liable_T_cnstr$p.value)
uncnstr_anno <- get_annotation(vB_m_v_m_liable_T_uncnstr$p.value)
dev.new(width=13,height=6,noRStudioGD = TRUE)
p3_3 <- ggplot(d_merged,aes(x=factor(scen_name),y=vB_m_v_m_liable, fill=agent_name)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(cnstr_anno[1]),unlist(uncnstr_anno[1])), textsize=7.5)
p3_3 <- p3_3 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Manufacturer Liability") +
  xlab ("Scenario Type") + ylab ("Measure") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=18, hjust=0.5)) +
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14), legend.position="top")+
  labs(fill='')+
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3_3

## (4) VA counterfactual
cnstr_anno <- get_annotation(vA_cntrfctl_T_cnstr$p.value)
uncnstr_anno <- get_annotation(vA_cntrfctl_T_uncnstr$p.value)
dev.new(width=13,height=6,noRStudioGD = TRUE)
p3_4 <- ggplot(d_merged,aes(x=factor(scen_name),y=vA_cntrfctl, fill=agent_name)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(cnstr_anno[1]),unlist(uncnstr_anno[1])), textsize=7.5)
p3_4 <- p3_4 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Consider Veh. A Counterfactual") +
  xlab ("Scenario Type") + ylab ("Measure") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=18, hjust=0.5)) +
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14), legend.position="top")+
  labs(fill='')+
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3_4

## (5) VB counterfactual
cnstr_anno <- get_annotation(vB_cntrfctl_T_cnstr$p.value)
uncnstr_anno <- get_annotation(vB_cntrfctl_T_uncnstr$p.value)
dev.new(width=13,height=6,noRStudioGD = TRUE)
p3_5 <- ggplot(d_merged,aes(x=factor(scen_name),y=vB_cntrfctl, fill=agent_name)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(cnstr_anno[1]),unlist(uncnstr_anno[1])), textsize=7.5)
p3_5 <- p3_5 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Consider Veh. B Counterfactual") +
  xlab ("Scenario Type") + ylab ("Measure") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=18, hjust=0.5)) +
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14), legend.position="top")+
  labs(fill='')+
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3_5

## (6) Capability to Avoid
cnstr_anno <- get_annotation(avoid_T_cnstr$p.value)
uncnstr_anno <- get_annotation(avoid_T_uncnstr$p.value)
dev.new(width=13,height=6,noRStudioGD = TRUE)
p3_6 <- ggplot(d_merged,aes(x=factor(scen_name),y=avoid, fill=agent_name)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(cnstr_anno[1]),unlist(uncnstr_anno[1])), textsize=7.5)
p3_6 <- p3_6 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Capability to Avoid") +
  xlab ("Scenario Type") + ylab ("Measure") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=18, hjust=0.5)) +
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14), legend.position="top")+
  labs(fill='')+
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p3_6

## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
t_names <- c("Constrained", "Unconstrained")

## (1) VA driver liable
annotations <- get_annotation(vA_liable_T_scen$p.value)
p2_1 <- ggplot(d_merged,aes(x=factor(scen_name),y=vA_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncnstr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p2_1 <- p2_1 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
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
p2_1

## (2) VB manufacturer/driver liability
annotations <- get_annotation(vB_m_v_d_liable_T_scen$p.value)
p2_2 <- ggplot(d_merged,aes(x=factor(scen_name),y=vB_m_v_d_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncnstr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p2_2 <- p2_2 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Manufacturer\nor Driver Liability") +
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
p2_2

## (3) VB manufacturer liability
annotations <- get_annotation(vB_m_v_m_liable_T_scen$p.value)
p2_3 <- ggplot(d_merged,aes(x=factor(scen_name),y=vB_m_v_m_liable)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncnstr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p2_3 <- p2_3 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
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
p2_3

## (4) VA counterfactual
annotations <- get_annotation(vA_cntrfctl_T_scen$p.value)
p2_4 <- ggplot(d_merged,aes(x=factor(scen_name),y=vA_cntrfctl)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncnstr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p2_4 <- p2_4 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
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
p2_4

## (5) VB Counterfactual
annotations <- get_annotation(vB_cntrfctl_T_scen$p.value)
p2_5 <- ggplot(d_merged,aes(x=factor(scen_name),y=vB_cntrfctl)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncnstr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p2_5 <- p2_5 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
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
p2_5

## (6) Capability to Avoid
annotations <- get_annotation(avoid_T_scen$p.value)
p2_6 <- ggplot(d_merged,aes(x=factor(scen_name),y=avoid)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c("cnstr", "uncnstr")), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p2_6 <- p2_6 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
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
p2_6

## PLOT SERIES 2
dev.new(width=10,height=8,noRStudioGD = TRUE)
figure2 <- ggarrange(p2_1, p2_2, p2_3, p2_4, p2_5, p2_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure2 <- annotate_figure(figure2,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Scenario Type", color="black", face ="plain",size=16)) 

plot(figure2)

write.csv(d_merged, 'd_spss.csv')

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================