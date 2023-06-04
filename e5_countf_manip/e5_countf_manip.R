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
               'Hmisc'
)

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
# if importing from Qualtrics: (i) export data as numeric values, and (ii) delete rows 2 and 3 of the .csv file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
d <- read.csv('e5_countf_manip_300_coded.csv')

## explore dataframe: 
dim <- dim(d) # will provide dimensions of the dataframe by row [1] and column [2]
colnames(d) # will provide all column names
summary(d)

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d) # number of participants should decrease after attention exclusions

## get number of participants BEFORE exclusions: 
n_original <- dim(d)[1] # extracting number of rows only, not columns

## perform comprehension exclusions separately for AV and HDV: 
# this will remove responses from the dataframe that failed comprehension checks (i.e., "2")
# countf_cat: 1 - human could have done more; 2 - no change or AV better
d <- subset(d, (d$comp1 == 1 & d$comp_accident == 1 & d$countf_cat != 0))
dim(d) # number of participants should decrease after comprehension exclusions

## get number of participants AFTER exclusions: 
n_final <- dim(d)[1] # extracting number of rows only, not columns
percent_excluded <- (n_original - n_final)/n_original

## get mean age and gender:
# mean(d$age, na.rm = TRUE) # removing NAs from the dataframe before computing mean 
mean_age = mean(as.numeric(d$age), na.rm = TRUE)
table(d$gender)[2]/sum(table(d$gender))

d <- d[,19:47]

d$vA_liable_AV_1 <- as.numeric(d$vA_liable_AV_1)
d$vB_mnfctr_liable_AV_1 <- as.numeric(d$vB_mnfctr_liable_AV_1)
d$avoid_AV_1 <- as.numeric(d$avoid_AV_1)

## ================================================================================================================
##                                             DATA ANALYSIS - T-TESTS               
## ================================================================================================================

## run t-tests
table(d$countf_cat)[1]/sum(table(d$countf_cat))
table(d$countf_cat)[2]/sum(table(d$countf_cat))
tapply(d$vB_mnfctr_liable_AV_1, d$countf_cat, mean)

## (1) LIABLE VEHICLE A DRIVER
vA_liable_T <- t.test(vA_liable_AV_1 ~ countf_cat, data = d, paired = FALSE) 
vA_liable_T
cohen.d(d$vA_liable_AV_1, d$countf_cat)

## (2) LIABLE VEHICLE B MANUFACTURER
vB_mnfctr_liable_T <- t.test(vB_mnfctr_liable_AV_1 ~ countf_cat, data = d, paired = FALSE) 
vB_mnfctr_liable_T
cohen.d(d$vB_mnfctr_liable_AV_1, d$countf_cat)

## (3) VEHICLE B CAN AVOID
avoid_T <- t.test(avoid_AV_1 ~ countf_cat, data = d, paired = FALSE) 
avoid_T
cohen.d(d$avoid_AV_1, d$countf_cat)

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

# SINGLE MEDIATION
process(data = d, y = "vB_mnfctr_liable_AV_1", x = "countf_cat",
        m ="avoid_AV_1", model = 4, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

## ================================================================================================================
##                                              PLOTTING MAIN FIGURES                 
## ================================================================================================================

## plotting all measures
t_names <- c("Yes", "No")
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
p1_1 <- ggplot(d,aes(x=factor(countf_cat),y=vA_liable_AV_1)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c(1, 2)), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

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

## (2) VB manufacturer liabile
annotations <- get_annotation(vB_mnfctr_liable_T$p.value)
p1_2 <- ggplot(d,aes(x=factor(countf_cat),y=vB_mnfctr_liable_AV_1)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c(1,2)), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_2 <- p1_2 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Veh. B Manufacturer\nLiability") +
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

## (3) Capability to Avoid
annotations <- get_annotation(avoid_T$p.value)
p1_3 <- ggplot(d,aes(x=factor(countf_cat),y=avoid_AV_1)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(comparisons = list(c(1,2)), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))

p1_3 <- p1_3 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Capability to Avoid") +
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
p1_3


## PLOT SERIES 1

dev.new(width=13,height=6,noRStudioGD = TRUE)
figure1 <- ggarrange(p1_1, p1_2, p1_3, nrow=1,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure1 <- annotate_figure(figure1,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                bottom = text_grob("Counterfactual Mentions Superior Human Driver", color="black", face ="plain",size=16))
plot(figure1)


#figure1 <- ggarrange(p1_2, p1_3, p1_5, nrow=1,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
#annotate_figure(figure1,left = text_grob("Mean Rating", color="black", face ="plain",size=16, rot=90),
#                bottom = text_grob("Scenario Condition", color="black", face ="plain",size=16)) 

plot(figure1)

write.csv(d_merged, 'd_spss.csv')

## ================================================================================================================
##                                                  Comparing E6 to E1-3                
## ================================================================================================================

#1=AV, 2=HDV

d$exp <- 6
names(d)[names(d)=="countf_cat"] <- "cond_n"
names(d)[names(d)=="vB_sue_AV_1"] <- "vB_sue"

#read in data
d1 <- read.csv('e1_processed.csv')
d1$exp <- 1
d2 <- read.csv('e2_processed.csv')
d2$exp <- 2
d3 <- read.csv('e3_processed.csv')
d3$exp <- 3

colnames(d[,c(12,14)])
colnames(d1[,c(16,4)])
colnames(d2[,c(17,4)])
colnames(d3[,c(17,4)])

#bind data
d_final <- rbind(d[,c(12,14,31)], d1[,c(16,4,17)], d2[,c(17,4,19)], d3[,c(17,4,18)])

d_f_s <- subset(d_final, d_final$cond_n == 1)

mod_6_1 <- t.test(d_f_s$vB_sue[d_f_s$exp==6], d_f_s$vB_sue[d_f_s$exp==1], paired = FALSE)  
mod_6_1

mod_6_2 <- t.test(d_f_s$vB_sue[d_f_s$exp==6], d_f_s$vB_sue[d_f_s$exp==2], paired = FALSE)  
mod_6_2

mod_6_3 <- t.test(d_f_s$vB_sue[d_f_s$exp==6], d_f_s$vB_sue[d_f_s$exp==3], paired = FALSE)  
mod_6_3

t_names <- c("Study 1", "Study 2", "Study 3", "Study 5")

p1_2_1 <- ggplot(d_f_s,aes(x=factor(exp),y=vB_sue)) +  
  theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  geom_signif(y_position = c(99,99), comparisons = list(c(1, 4)), annotation="***", textsize = 6)+
  geom_signif(y_position = c(108,106), comparisons = list(c(2, 4)), annotation="***", textsize = 6)+
  geom_signif(y_position = c(119,116), comparisons = list(c(3, 4)), annotation="***", textsize = 6)+
  coord_cartesian(ylim=c(1,125)) 

p1_2_1 <- p1_2_1 + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=t_names) +
  ggtitle("Sue Veh. B Manufacturer") +
  
  xlab ("Experiment") + ylab ("Mean Agreement") +
  theme_classic() +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.title = element_text(size=18)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(plot.title = element_text(size=20, hjust=0.5)) +
  geom_violin(width=0.9, alpha=0.38, size=0.75) +  
  geom_sina(alpha=0.6, size=0.95, color = "#999999") +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               size=0.4, 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_cl_boot", color = "black", 
               position = position_dodge(width = 0.9),
               geom="errorbar", width = 0.2)
p1_2_1

dev.new(width=13,height=6,noRStudioGD = TRUE)
p1_2_1

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================