## clear workspace
rm(list = ls()) 

options(download.file.method="libcurl") # sets method for downloading files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to current directory

source("../common.R") # install packages; import common plotting functions

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
d_raw <- read.csv('e2_counterf_gen.csv')

## explore dataframe: 
dim(d_raw) 

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d_raw <- subset(d_raw, (d_raw$att1 == 2 & d_raw$att2 == 2))
dim(d_raw) 

## new columns to label agent (av/hdv), scenario (cnstr/uncnstr), and condition (combination of the 2)
d_raw$agent_cond <- ifelse(is.na(d_raw$vA_sue_AV_1), "hdv", "av")
d_raw$cond <- paste(d_raw$agent, d_raw$countf, sep="_")

## get number of participants before comp exclusions
n_orig_all <- dim(d_raw)[1]
n_orig <- as.list(table(d_raw$cond))

## perform comp exclusions
d_raw <- subset(d_raw, comp_accident == 1) # include only comp check 2 passes
d_raw <- subset(d_raw, cruise_familiarity == 2)
d_raw <- subset(d_raw, (agent_cond == "av" & comp1 == 1) | (agent_cond == "hdv" & comp1 == 2))

#coder reliability
cor.test(d_raw$coder_1, d_raw$coder_2)
d_raw <- subset(d_raw, (coder_1 == coder_2))

cor.test(d_raw$category_coder1[!is.na(d_raw$category_coder1)], 
         d_raw$category_coder2[!is.na(d_raw$category_coder2)])

## get number of participants AFTER exclusions: 
n_ss <- dim(d_raw)[1]
percent_excl_all <- (n_orig_all - n_ss)/n_orig_all; percent_excl_all
n_excl_all <- n_orig_all - n_ss; n_excl_all
n_excl_all/n_orig_all
table(d_raw$cond)

## get mean age and gender:
mean_age = mean(as.numeric(d_raw$age), na.rm = TRUE); mean_age
gender_f = table(d_raw$gender)["2"]/sum(table(d_raw$gender)); gender_f

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_merged <- array(dim=c(0, 7))
colnames(d_merged) <- c('vA_sue', 'vB_mnfctr_liab', 'vB_mnfctr_liab_2', 'vA_cntrfctl', 'vB_cntrfctl_1', 'vB_cntrfctl_2',
                        'counterfactual_content')

d_merged <- as.data.frame(d_merged, stringsAsFactors=FALSE) 

## select only the used columns
cols_subset <- c(26:31, 36:41)
for(i in 1:n_ss) {
    curr <- d_raw[i,cols_subset][!is.na(d_raw[i,cols_subset])] #for a given row, get only the non NA values
    d_merged[i,1:6] <- as.numeric(curr[curr!= ""]) 
    d_merged[i,7] <- ifelse(all(d_raw[i,42:43] == ""), d_raw[i,42:43][1], d_raw[i,42:43][d_raw[i,42:43] != ""])
}

# make columns numeric
d_merged[,1:6] <- lapply(d_merged[,1:6], as.numeric)
d <- cbind(d_merged, d_raw[,44:ncol(d_raw)])

# agent_n where av=1, human=2; intv_n where yes=1, no=2
d$agent_n <- ifelse(d$agent_cond=="av", 1, 2)

# reliability
cor(d$vB_mnfctr_liab, d$vB_mnfctr_liab_2)
d$vB_liability <-  (d$vB_mnfctr_liab + d$vB_mnfctr_liab_2)/2 ###

cor(d$vB_cntrfctl_1, d$vB_cntrfctl_2)
d$counterf_med <- (d$vB_cntrfctl_1 + d$vB_cntrfctl_2)/2 ###

#discriminant validity
library(lavaan)
library(semTools)

countf.model <- ' liability   =~ vB_mnfctr_liab + vB_mnfctr_liab_2
                  counterfactual  =~ vB_cntrfctl_1 + vB_cntrfctl_2 '
 
htmt(countf.model, d)

## save covariance matrix
countf.cov <- cov(d[, c(2:3, 5:6)])

## HTMT using arithmetic mean
htmt(countf.model, sample.cov = countf.cov, htmt2 = FALSE)

## ================================================================================================================
##                                              PLOTTING 2X2 FIGURE                 
## ================================================================================================================

fill_labels <- c("AV", "HDV")
sig_comparisons <- c("av", "hdv")
d$code <- d$coder_1

#basic liability effect
p_val = t.test(vB_liability ~ agent_cond, data = d)$p.value
p0_1 <- plot_std(d, x=agent_cond, y=vB_liability, p_val, 
                 title="Firm Liability", fill_labels, sig_comparisons)

#basic counterfactual effect
p_val = t.test(counterf_med ~ agent_cond, data = d)$p.value
p0_2 <- plot_std(d, x=agent_cond, y=counterf_med, p_val, 
                 title="Counterfactual Relevance", fill_labels, sig_comparisons)

## Understanding why counterfactual means are flat
p_val_L = t.test(counterf_med ~ code, data = d[d$agent_cond == "av", ])$p.value
p_val_R = t.test(counterf_med ~ code, data = d[d$agent_cond == "hdv", ])$p.value

code_labels <- c("1", "2")

p0_3 <- plot_2x2(d, x=agent_cond, y=counterf_med, fill=as.factor(code), p_val_L, p_val_R, 
                 title="Counterfactual Relevance", c("AV", "HDV"), code_labels)
p0_3

figure0 <- ggarrange(p0_1, p0_2, p0_3, nrow=1,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure0 <- annotate_figure(figure0,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 
figure0

#---

## 2x2 liability based on coutnerfactual codes
p_val_L = t.test(vB_liability ~ agent_cond, data = d[d$code == 1, ])$p.value
p_val_R = t.test(vB_liability ~ agent_cond, data = d[d$code == 2, ])$p.value

p1_1 <- plot_2x2(d, x=code, y=vB_liability, fill=agent_cond, p_val_L, p_val_R, 
                 title="Firm Liability", c("No", "Yes"), fill_labels)
p1_1


## 2x2 liability based on coutnerfactual codes
p_val_L = t.test(counterf_med ~ agent_cond, data = d[d$code == 1, ])$p.value
p_val_R = t.test(counterf_med ~ agent_cond, data = d[d$code == 2, ])$p.value

p1_2 <- plot_2x2(d, x=code, y=counterf_med, fill=agent_cond, p_val_L, p_val_R, 
                 title="Counterfactual Relevance", c("No", "Yes"), fill_labels)
p1_2

figure2 <- ggarrange(p1_1, p1_2, nrow=1,ncol=2,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure2 <- annotate_figure(figure2,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Counterfactual Mentions Not-At-Fault Vehicle", color="black", face ="plain",size=16)) 

#png(file = "../plots/e2_plot.png", width = 2*900, height = 2*700, res = 200)  # width and height are in inches
plot(figure2)
#dev.off()

## ================================================================================================================
##                                         DATA ANALYSIS - ANOVA              
## ================================================================================================================

## Proportion of mentions of not-at-fault vehicle
table(d$agent_cond, d$code)
N <- sum(table(d$agent_cond, d$code))
prop.table(table(d$agent_cond, d$code), margin = 1) * 100
chisq.test(as.factor(d$code), d$agent_cond)


## Veh A Sue 
t.test(vA_sue ~ agent_cond, data = d, subset = code == 2)
cohen.d(d$vA_sue[d$code==2], d$agent_cond[d$code==2])

t.test(vA_sue ~ agent_cond, data = d, subset = code == 1)
cohen.d(d$vA_sue[d$code==1], d$agent_cond[d$code==1])


## Veh A Counterfactual 
t.test(vA_cntrfctl ~ agent_cond, data = d, subset = code == 2)
cohen.d(d$vA_cntrfctl[d$code==2], d$agent_cond[d$code==2])

t.test(vA_cntrfctl ~ agent_cond, data = d, subset = code == 1)
cohen.d(d$vA_cntrfctl[d$code==1], d$agent_cond[d$code==1])


## Veh B Liability
liab_coded_mod <- aov(vB_liability ~ as.factor(agent_n) * as.factor(code), data = d)
summary(liab_coded_mod)
anova_stats(liab_coded_mod)

aggregate(vB_liability ~ agent_cond, data = d, FUN = mean)
aggregate(vB_liability ~ code, data = d, FUN = mean)

#t-tests
t.test(vB_liability ~ agent_cond, data = d, subset = code == 2)
cohen.d(d$vB_liability[d$code==2], d$agent_cond[d$code==2])

t.test(vB_liability ~ agent_cond, data = d, subset = code == 1)
cohen.d(d$vB_liability[d$code==1], d$agent_cond[d$code==1])


## Veh B Counterfactual rel
countf_coded_mod <- aov(counterf_med ~ as.factor(agent_n) * as.factor(code), data = d)
summary(countf_coded_mod)
anova_stats(countf_coded_mod)

aggregate(counterf_med ~ agent_cond, data = d, FUN = mean)
aggregate(counterf_med ~ code, data = d, FUN = mean)

#t-tests
t.test(counterf_med ~ code, data = d, subset = agent_cond == "hdv")
cohen.d(d$counterf_med[d$agent_cond=="hdv"], d$code[d$agent_cond=="hdv"])

t.test(counterf_med ~ code, data = d, subset = agent_cond == "av")
cohen.d(d$counterf_med[d$agent_cond=="av"], d$code[d$agent_cond=="av"])


## Coding of exploratory analysis
table(d$category_coder1)
prop.table(table(d$agent_cond, d$category_coder1), margin = 1) * 100
observed_counts <- c(23,40,13)
expected_counts <- rep(N/3, length(observed_counts))
chisq.test(observed_counts, p = rep(1/3, length(observed_counts)))

aggregate(vB_liability ~ category_coder1, data = d, FUN = mean)
liab_coded_mod <- aov(vB_liability ~ as.factor(category_coder1), data = d)
summary(liab_coded_mod)
anova_stats(liab_coded_mod)

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

mediation <- FALSE #change to true if you want to run this code

if(mediation) {
  source("../process.R")
  
  # test age as moderator
  summary(lm(vB_liability ~ agent_n*age, data=d))
  
  # Simple mediation
  process(data = d, y = "vB_liability", x = "agent_n",
          m =c("counterf_med"), model = 4, effsize =1, total =1, stand =1,
          contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
  
  # MODERATED MEDIATION
  # the effect of intervention on A path (7)
  process(data = d, y = "vB_liability", x = "agent_n",
          m =c("counterf_med"), w = "code", model = 7, effsize =1, total =1, stand =1,
          contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
}

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================