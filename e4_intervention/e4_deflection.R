## clear workspace
rm(list = ls()) 

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
               'Hmisc',
               'sjstats'
)

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
# if importing from Qualtrics: (i) export data as numeric values, and (ii) delete rows 2 and 3 of the .csv file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
d <- read.csv('e4_deflection_1350.csv')

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
d_AV_intv <- subset(d, (d$FL_92_DO == "FL_94"))
d_HDV_intv <- subset(d, (d$FL_98_DO == "FL_100"))
d_AV_no_intv <- subset(d, (d$FL_92_DO == "FL_93"))
d_HDV_no_intv <- subset(d, (d$FL_98_DO == "FL_99"))

## get number of participants BEFORE exclusions: 
n_original <- dim(d)[1] # extracting number of rows only, not columns
n_original_AV_intv <- dim(d_AV_intv)[1]
n_original_HDV_intv <- dim(d_HDV_intv)[1]
n_original_AV_no_intv <- dim(d_AV_no_intv)[1]
n_original_HDV_no_intv <- dim(d_HDV_no_intv)[1]

## perform comprehension exclusions separately for AV and HDV: 
# this will remove responses from the dataframe that failed comprehension checks (i.e., "2")
d_AV_intv <- subset(d_AV_intv, (d_AV_intv$comp1 == 1 & d_AV_intv$comp_accident == 1))
d_HDV_intv <- subset(d_HDV_intv, (d_HDV_intv$comp1 == 2 & d_HDV_intv$comp_accident == 1))
d_AV_no_intv <- subset(d_AV_no_intv, (d_AV_no_intv$comp1 == 1 & d_AV_no_intv$comp_accident == 1))
d_HDV_no_intv <- subset(d_HDV_no_intv, (d_HDV_no_intv$comp1 == 2 & d_HDV_no_intv$comp_accident == 1))

## perform exclusions for intervention
d_AV_intv <- subset(d_AV_intv, (d_AV_intv$intv_check_AV == "1,2"))
d_HDV_intv <- subset(d_HDV_intv, (d_HDV_intv$intv_check == "1,2"))

## get number of participants AFTER exclusions: 
n_final_AV_intv <- dim(d_AV_intv)[1] # extracting number of rows only, not columns
n_final_HDV_intv <- dim(d_HDV_intv)[1]
n_final_AV_no_intv <- dim(d_AV_no_intv)[1]
n_final_HDV_no_intv <- dim(d_HDV_no_intv)[1]
n_final <- n_final_AV_intv + n_final_HDV_intv + n_final_AV_no_intv + n_final_HDV_no_intv
percent_excluded <- (n_original - n_final)/n_original
percent_excluded_AV_intv <- (n_original_AV_intv - n_final_AV_intv)/n_original_AV_intv
percent_excluded_HDV_intv <- (n_original_HDV_intv - n_final_HDV_intv)/n_original_HDV_intv
percent_excluded_AV_no_intv <- (n_original_AV_no_intv - n_final_AV_no_intv)/n_original_AV_no_intv
percent_excluded_HDV_no_intv <- (n_original_HDV_no_intv - n_final_HDV_no_intv)/n_original_HDV_no_intv

## remove unused columns (other condition and click info) according to condition
d_AV_intv <- d_AV_intv[-c(21:33,34:42,48:76)] # first=prior columns, second=click, third=all other columns
d_HDV_intv <- d_HDV_intv[-c(21:61,62:70)] # first=prior columns, second=click info
d_AV_no_intv <- d_AV_no_intv[-c(21:28,34:76)] # first=click columns, second=later columns
d_HDV_no_intv <- d_HDV_no_intv[-c(21:47,48:55,62:76)] # first=prior columns, second=click info

## duplicate AV condition vB liability column to match with HDV condition driver liability col
d_AV_intv$vB_mf_sue_AV_intv_2 <- d_AV_intv$vB_mf_sue_AV_intv_1
d_AV_intv <- d_AV_intv %>% relocate(vB_mf_sue_AV_intv_2, .after=vB_mf_sue_AV_intv_1)
d_AV_no_intv$vB_mf_sue_AV_2 <- d_AV_no_intv$vB_mf_sue_AV_1
d_AV_no_intv <- d_AV_no_intv %>% relocate(vB_mf_sue_AV_2, .after=vB_mf_sue_AV_1)

## get mean age and gender:
mean_age = mean(as.numeric(d$age), na.rm = TRUE) # removing NAs from the dataframe before computing mean 
gender = table(d$gender)["1"]/sum(table(d$gender)) # percent male

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_subset <- array(dim=c(n_final, 11))
colnames(d_subset) <- c('agent_name', 'intv_appld', 'vA_sue', 'vB_m_v_d_sue', 'vB_m_v_m_sue',
                        'vA_cntrfctl', 'vB_cntrfctl', 'avoid', 'comp1', 'comp2', 'age')
d_subset <- as.data.frame(d_subset, stringsAsFactors=FALSE) 

## extract good data from the middle part of raw data in AV + intervention
for(i in 1:n_final_AV_intv) {
  curr <- d_AV_intv[i,21:28][!is.na(d_AV_intv[i,21:28])] # for a given row, get only the non-NA values
  d_subset[i,3:10] <- as.numeric(curr[curr!= ""]) # and only the non-empty values
  #d_subset[i,1] <- d_AV[i,43][!is.na(d_AV[i,43])]
  d_subset[i,11] <- as.numeric(d_AV_intv$age[i])
  d_subset[i,1] <- "av"
  d_subset[i,2] <- "yes"
}

## extract good data from the middle part of raw data in HDV + intervention
for(i in 1:n_final_HDV_intv) {
  j = i+n_final_AV_intv
  curr <- d_HDV_intv[i,21:28][!is.na(d_HDV_intv[i,21:28])] # for a given row, get only the non-NA values
  d_subset[j,3:10] <- as.numeric(curr) # and only the non-empty values
  d_subset[j,11] <- as.numeric(d_HDV_intv$age[i])
  d_subset[j,1] <- "human"
  d_subset[j,2] <- "yes"
}

## extract good data from the middle part of raw data in AV + no intervention
for(i in 1:n_final_AV_no_intv) {
  j = i+n_final_AV_intv+n_final_HDV_intv
  curr <- d_AV_no_intv[i,21:28][!is.na(d_AV_no_intv[i,21:28])] # for a given row, get only the non-NA values
  d_subset[j,3:10] <- as.numeric(curr) # and only the non-empty values
  d_subset[j,11] <- as.numeric(d_AV_no_intv$age[i])
  d_subset[j,1] <- "av"
  d_subset[j,2] <- "no"
}

## extract good data from the middle part of raw data in HDV + no intervention
for(i in 1:n_final_HDV_no_intv) {
  j = i+n_final_AV_intv+n_final_HDV_intv+n_final_AV_no_intv
  curr <- d_HDV_no_intv[i,21:28][!is.na(d_HDV_no_intv[i,21:28])] # for a given row, get only the non-NA values
  d_subset[j,3:10] <- as.numeric(curr) # and only the non-empty values
  d_subset[j,11] <- as.numeric(d_HDV_no_intv$age[i])
  d_subset[j,1] <- "human"
  d_subset[j,2] <- "no"
}

## just to keep the df names straight for next section
d_merged <- d_subset

# agent_n where av=1, human=2; intv_n where yes=1, no=2
d_merged$agent_n <- ifelse(d_merged$agent_name=="av", 1, 2)
d_merged$intv_n <- ifelse(d_merged$intv_appld=="yes", 1, 2)

## ================================================================================================================
##                                      DATA ANALYSIS - ANOVA & T-TESTS               
## ================================================================================================================

#### get summary statistics
d_merged %>%
  group_by(intv_appld, agent_name) %>%
  get_summary_stats(vB_m_v_d_sue, type = "mean_sd")
d_merged %>%
  group_by(intv_appld, agent_name) %>%
  get_summary_stats(vB_m_v_m_sue, type = "mean_sd")
d_merged %>%
  group_by(intv_appld, agent_name) %>%
  get_summary_stats(vB_cntrfctl, type = "mean_sd")
d_merged %>%
  group_by(intv_appld, agent_name) %>%
  get_summary_stats(avoid, type = "mean_sd")

##### ----- ANOVA -----
### manufacturer vs driver
m_v_d_mod <- aov(vB_m_v_d_sue ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(m_v_d_mod)
anova_stats(m_v_d_mod)

### manufacturer vs manufacturer
m_v_m_mod <- aov(vB_m_v_m_sue ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(m_v_m_mod)
anova_stats(m_v_m_mod)

### counterfactual
countf_mod <- aov(vB_cntrfctl ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(countf_mod)
anova_stats(countf_mod)

### avoid
avoid_mod <- aov(avoid ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(avoid_mod)
anova_stats(avoid_mod)

##### ----- T-TESTS -----

### (1) sue vehicle A driver
# does measure depend on agent
vA_sue_T_agent <- t.test(vA_sue ~ agent_name, data = d_merged, paired = FALSE) 

# does measure depend on intervention
vA_sue_T_scen <- t.test(vA_sue ~ intv_appld, data = d_merged, paired = FALSE) 

# when intervention, do av and human differ
vA_sue_T_intv <- t.test(d_merged$vA_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                            d_merged$vA_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vA_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
       d_merged$vA_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"])
# when no intervention, do av and human differ
vA_sue_T_no_intv <- t.test(d_merged$vA_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
                              d_merged$vA_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vA_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
       d_merged$vA_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"])


### (2) sue vehicle B manufacturer VS driver
# does measure depend on agent
vB_m_v_d_sue_T_agent <- t.test(vB_m_v_d_sue ~ agent_name, data = d_merged, paired = FALSE) 

# does measure depend on intervention
vB_m_v_d_sue_T_scen <- t.test(vB_m_v_d_sue ~ intv_appld, data = d_merged, paired = FALSE) 

# when intervention, does av and human differ
vB_m_v_d_sue_T_intv <- t.test(d_merged$vB_m_v_d_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                                  d_merged$vB_m_v_d_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], paired=FALSE)

cohen.d(d_merged$vB_m_v_d_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
       d_merged$vB_m_v_d_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"])

# when no intervention, does av and human differ
vB_m_v_d_sue_T_no_intv <- t.test(d_merged$vB_m_v_d_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
                                    d_merged$vB_m_v_d_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)

cohen.d(d_merged[d_merged$intv_appld=="no", ]$vB_m_v_d_sue, d_merged[d_merged$intv_appld=="no", ]$agent_name)
cohen.d(d_merged$vB_m_v_d_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
       d_merged$vB_m_v_d_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"])

# when av, do intervention differ
vB_m_v_d_sue_T_av <- t.test(d_merged$vB_m_v_d_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                      d_merged$vB_m_v_d_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], paired=FALSE)

# when hdv, do intervention differ
vB_m_v_d_sue_T_hdv <- t.test(d_merged$vB_m_v_d_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], 
                       d_merged$vB_m_v_d_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)


### (3) sue vehicle B manufacturer VS manufacturer
# t-test; does measure depend on agent
vB_m_v_m_sue_T_agent <- t.test(vB_m_v_m_sue ~ agent_name, data = d_merged, paired = FALSE) 

# t-test; does measure depend on intervention
vB_m_v_m_sue_T_scen <- t.test(vB_m_v_m_sue ~ intv_appld, data = d_merged, paired = FALSE) 

# when intervention, do av and human differ
vB_m_v_m_sue_T_intv <- t.test(d_merged$vB_m_v_m_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                                  d_merged$vB_m_v_m_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vB_m_v_m_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
       d_merged$vB_m_v_m_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"])

cohen.d(d_merged[d_merged$intv_appld=="yes", ]$vB_m_v_m_sue, d_merged[d_merged$intv_appld=="yes", ]$agent_name)

# when no intervention, do av and human differ
vB_m_v_m_sue_T_no_intv <- t.test(d_merged$vB_m_v_m_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
                                    d_merged$vB_m_v_m_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)

cohen.d(d_merged$vB_m_v_m_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
       d_merged$vB_m_v_m_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"])

# when av, do intervention differ
vB_m_v_m_sue_T_av <- t.test(d_merged$vB_m_v_m_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                            d_merged$vB_m_v_m_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], paired=FALSE)

# when hdv, do intervention differ
vB_m_v_m_sue_T_hdv <- t.test(d_merged$vB_m_v_m_sue[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], 
                             d_merged$vB_m_v_m_sue[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)

### (4) vehicle A counterfactual
# does measure depend on agent
vA_cf_T_agent <- t.test(vA_cntrfctl ~ agent_name, data = d_merged, paired = FALSE) 

# does measure depend on intervention
vA_cf_T_scen <- t.test(vA_cntrfctl ~ intv_appld, data = d_merged, paired = FALSE) 

# when intervention, do av and human differ
vA_cf_T_intv <- t.test(d_merged$vA_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                       d_merged$vA_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vA_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
       d_merged$vA_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"])

# when no intervention, do av and human differ
vA_cf_T_no_intv <- t.test(d_merged$vA_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
                          d_merged$vA_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vA_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
       d_merged$vA_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "human"])


### (5) vehicle B counterfactual
# does measure depend on agent
vB_cf_T_agent <- t.test(vB_cntrfctl ~ agent_name, data = d_merged, paired = FALSE) 

# does measure depend on intervention
vB_cf_T_scen <- t.test(vB_cntrfctl ~ intv_appld, data = d_merged, paired = FALSE) 

# when intervention, do av and human differ
vB_cf_T_intv <- t.test(d_merged$vB_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                       d_merged$vB_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vB_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
       d_merged$vB_cntrfctl[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"])

# when no intervention, do av and human differ
vB_cf_T_no_intv <- t.test(d_merged$vB_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
                          d_merged$vB_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$vB_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
      d_merged$vB_cntrfctl[d_merged$intv_appld=="no" & d_merged$agent_name == "human"])


### (6) vehicle B can avoid
# does measure depend on agent
avoid_T_agent <- t.test(avoid ~ agent_name, data = d_merged, paired = FALSE) 

# does measure depend on intervention
avoid_T_scen <- t.test(avoid ~ intv_appld, data = d_merged, paired = FALSE) 

# when intervention, do av and human differ
avoid_T_intv <- t.test(d_merged$avoid[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
                        d_merged$avoid[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$avoid[d_merged$intv_appld=="yes" & d_merged$agent_name == "av"], 
       d_merged$avoid[d_merged$intv_appld=="yes" & d_merged$agent_name == "human"])

# when no intervention, do av and human differ
avoid_T_no_intv <- t.test(d_merged$avoid[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
                        d_merged$avoid[d_merged$intv_appld=="no" & d_merged$agent_name == "human"], paired=FALSE)
cohen.d(d_merged$avoid[d_merged$intv_appld=="no" & d_merged$agent_name == "av"], 
       d_merged$avoid[d_merged$intv_appld=="no" & d_merged$agent_name == "human"])


cor(d_merged[,3:8])

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

source("../process.R")

summary(lm(vB_m_v_m_sue ~ agent_n*age, data=d_merged))

# SERIAL MEDIATION
#  agent -> countfactual -> avoid -> sue
process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n",
        m =c("vB_cntrfctl", "avoid"), model = 6, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
process(data = d_merged, y = "vB_m_v_d_sue", x = "agent_n",
        m =c("vB_cntrfctl", "avoid"), model = 6, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

# MODERATED SERIAL MEDIATION
#  agent -> counterfactual -> avoid -> sue (moderated by intervention)
# 87 = B path, 83 = A path, 91 = center path
process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n",
        m =c("vB_cntrfctl", "avoid"), w = "intv_n", model = 83, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
process(data = d_merged, y = "vB_m_v_d_sue", x = "agent_n",
        m =c("vB_cntrfctl", "avoid"), w = "intv_n", model = 83, effsize =1, total =1, stand =1,
        contrast =1, boot = 10000 , modelbt = 1, seed = 654321)

## ================================================================================================================
##                                           PLOTTING MAIN FIGURE - BY AGENT                 
## ================================================================================================================

## x-axis = agent
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

plotting <- function (cond, x, y, title,fill_labels) {
  annotations <- get_annotation(cond$p.value)
  p <- ggplot(d_merged,aes(x=factor({{x}}),y={{y}})) +  
    theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
    geom_signif(comparisons = list(c(fill_labels)), annotation=unlist(annotations[1]), textsize = unlist(annotations[2]))
  p <- p + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_x_discrete(labels=t_names) +
    ggtitle(title) +
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
  return(p)
}
labels <- c("av", "human")

### (1) sue vehicle A driver
p1_1 <- plotting(vA_sue_T_agent,factor(agent_name),vA_sue,"Veh. A Driver Sue",labels)

## (2) sue vehicle B manufacturer VS driver
p1_2 <- plotting(vB_m_v_d_sue_T_agent,factor(agent_name),vB_m_v_d_sue,"Veh. B Manufacturer\nor Driver Sue",labels)

## (3) sue vehicle B manufacturer VS manufacturer
p1_3<-plotting(vB_m_v_m_sue_T_agent,factor(agent_name),vB_m_v_m_sue,"Veh. B Manufacturer Sue",labels)

## (4) vehicle A counterfactual
p1_4 <- plotting(vA_cf_T_agent,factor(agent_name),vA_cntrfctl,"Consider Veh. A Counterfactual",labels)

## (5) vehicle B counterfactual
p1_5 <- plotting(vB_cf_T_agent,factor(agent_name),vB_cntrfctl,"Consider Veh. B Counterfactual",labels)

## (6) vehicle B - capability to avoid
p1_6 <- plotting(avoid_T_agent,factor(agent_name),avoid,"Capability to Avoid",labels)

## PLOT SERIES 1 - BY AGENT
figure1 <- ggarrange(p1_1, p1_2, p1_3, p1_4, p1_5, p1_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure1 <- annotate_figure(figure1,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 

plot(figure1)

## ================================================================================================================
##                                  PLOTTING MAIN FIGURE - BY INTERVENTION                 
## ================================================================================================================

### x-axis = intervention
t_names <- c("No Intervention", "Intervention")
labels <- c("no","yes")

## (1) sue vehicle A driver
p2_1 <- plotting(vA_sue_T_scen,intv_appld,vA_sue,"Veh. A Driver Sue",labels)

## (2) sue vehicle B manufacturer VS driver
p2_2 <- plotting(vB_m_v_d_sue_T_scen,intv_appld,vB_m_v_d_sue,"Veh. B Manufacturer\nor Driver Sue",labels)

## (3) sue vehicle B manufacturer VS manufacturer
p2_3 <- plotting(vB_m_v_m_sue_T_scen,intv_appld,vB_m_v_m_sue,"Veh. B Manufacturer Sue",labels)

## (4) vehicle A counterfactual
p2_4 <- plotting(vA_cf_T_scen,intv_appld,vA_cntrfctl,"Consider Veh. A Counterfactual",labels)

## (5) vehicle B counterfactual
p2_5 <- plotting(vB_cf_T_scen,intv_appld,vB_cntrfctl,"Consider Veh. B Counterfactual",labels)

## (6) vehicle B - capability to avoid
p2_6 <- plotting(avoid_T_scen,intv_appld,avoid,"Capability to Avoid",labels)

## PLOT SERIES 2 - BY INTERVENTION
#dev.new(width=10,height=8,noRStudioGD = TRUE)
figure2 <- ggarrange(p2_1, p2_2, p2_3, p2_4, p2_5, p2_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure2 <- annotate_figure(figure2,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Intervention Applied", color="black", face ="plain",size=16)) 

plot(figure2)

## ================================================================================================================
##                                PLOTTING MAIN FIGURE - BY INTERVENTION --> AGENT                   
## ================================================================================================================

## x-axis = intervention
t_names <- c("No Intervention", "Intervention")
plotting_intervention <- function(cond_intv,cond_no_intv,y,title,fill_labels)
{
  intv_anno <- get_annotation(cond_intv$p.value)
  no_intv_anno <- get_annotation(cond_no_intv$p.value)
  p <- ggplot(d_merged,aes(x=factor(intv_appld),y={{y}}, fill=agent_name)) +  
    theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
    geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(intv_anno[1]),unlist(no_intv_anno[1])), textsize=6)
  p <- p + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_x_discrete(labels=t_names) +
    ggtitle(title) +
    xlab ("") + ylab ("") +
    scale_fill_discrete(labels=fill_labels) +
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
  return(p)
}
labels=c('AV', 'HDV')
# (1) sue vehicle A driver
p3_1<-plotting_intervention(vA_sue_T_intv,vA_sue_T_no_intv,vA_sue,"Veh. A Driver Sue",labels)
  
## (2) sue vehicle B manufacturer VS driver
p3_2<-plotting_intervention(vB_m_v_d_sue_T_intv,vB_m_v_d_sue_T_no_intv,vB_m_v_d_sue,"Veh. B Manufacturer\nor Driver Sue",labels)

## (3) sue vehicle B manufacturer VS manufacturer
p3_3 <- plotting_intervention(vB_m_v_m_sue_T_intv,vB_m_v_m_sue_T_no_intv,vB_m_v_m_sue,"Veh. B Manufacturer Sue",labels)

## (4) vA counterfactual
p3_4 <- plotting_intervention(vA_cf_T_intv,vA_cf_T_no_intv,vA_cntrfctl,"Consider Veh. A Counterfactual",labels)

## (5) vB counterfactual
p3_5 <- plotting_intervention(vB_cf_T_intv,vB_cf_T_no_intv,vB_cntrfctl,"Consider Veh. B Counterfactual",labels)

## (6) vehicle B - capability to avoid
p3_6 <- plotting_intervention(avoid_T_intv,avoid_T_no_intv,avoid,"Capability to Avoid",labels)

## PLOT SERIES 3 AND 4 - BY INTERVENTION --> AGENT
#dev.new(width=10,height=8,noRStudioGD = TRUE)
figure3 <- ggarrange(p3_1, p3_2, p3_3, p3_4, p3_5, p3_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure3 <- annotate_figure(figure3,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Intervention Applied", color="black", face ="plain",size=16)) 
plot(figure3)

#dev.new(width=10,height=5,noRStudioGD = TRUE)
figure4 <- ggarrange(p3_2, p3_3, nrow=1,ncol=2,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure4 <- annotate_figure(figure4,left = text_grob("Mean Agreement", color="black", face ="plain",size=18, rot=90),
                           bottom = text_grob("Intervention Applied", color="black", face ="plain",size=18)) 
plot(figure4)


## ================================================================================================================
##                                PLOTTING MAIN FIGURE - BY AGENT --> INTERVENTION                   
## ================================================================================================================

## x-axis = intervention
t_names <- c("No Intervention", "Intervention")

# (1) sue vehicle A driver
p4_1 <- plotting_intervention(vA_sue_T_av,vA_sue_T_hdv,vA_sue,"Veh. A Driver Sue",labels)

## (2) sue vehicle B manufacturer VS driver
p4_2 <- plotting_intervention(vB_m_v_d_sue_T_av,vB_m_v_d_sue_T_hdv,vB_m_v_d_sue,"Veh. B Manufacturer\nor Driver Sue",labels)

## (3) sue vehicle B manufacturer VS manufacturer
p4_3 <- plotting_intervention(vB_m_v_m_sue_T_av,vB_m_v_m_sue_T_hdv,vB_m_v_m_sue,"Veh. B Manufacturer Sue",labels)

## (4) vA counterfactual
p4_4 <- plotting_intervention(vA_cf_T_av,vA_cf_T_hdv,vA_cntrfctl,"Consider Veh. A Counterfactual",labels)

## (5) vB counterfactual
p4_5 <- plotting_intervention(vB_cf_T_av,vB_cf_T_hdv,vB_cntrfctl,"Consider Veh. B Counterfactual",labels)

## (6) vehicle B - capability to avoid
p4_6 <- plotting_intervention(avoid_T_av,avoid_T_hdv,avoid,"Capability to Avoid",labels)

## PLOT SERIES 3 AND 4 - BY INTERVENTION --> AGENT
#dev.new(width=10,height=8,noRStudioGD = TRUE)
figure5 <- ggarrange(p4_1, p4_2, p4_3, p4_4, p4_5, p4_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure5 <- annotate_figure(figure5, left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 
plot(figure5)

#dev.new(width=10,height=5,noRStudioGD = TRUE)
figure6 <- ggarrange(p4_2, p4_3, nrow=1,ncol=2,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure6 <- annotate_figure(figure6, left = text_grob("Mean Agreement", color="black", face ="plain",size=18, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=18)) 
plot(figure6)


## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================