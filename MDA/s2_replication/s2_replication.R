## clear workspace
rm(list = ls()) 

options(download.file.method="libcurl") # sets method for downloading files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to current directory

source("../../common.R") # install packages; import common plotting functions

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
d <- read.csv('s2_replication.csv')
dim(d) 

## perform attention exclusions: 
# this will remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d) 

## new column to label agent (av/hdv) condition
d$cond <- ifelse(d$FL_4_DO == "FL_39", "av", "hdv") #if FL_40 -> HDV

## get number of participants before comp exclusions
n_orig_all <- dim(d)[1]
n_orig <- as.list(table(d$cond))

## perform comp exclusions
d_clean <- d
d_clean <- subset(d_clean, comp_accident == 1) # include only comp check 2 passes
d_clean <- subset(d_clean, !(cond == "av" & comp1 != 1)) # remove comp check 1 fails for av
d_clean <- subset(d_clean, !(cond == "hdv" & comp1 != "")) # remove comp check 1 fails for hdv

## get number of participants AFTER exclusions: 
n_final_all <- dim(d_clean)[1]; n_final_all
percent_excl_all <- (n_orig_all - n_final_all)/n_orig_all
n_excl_all <- n_orig_all - n_final_all
n_final <- as.list(table(d_clean$cond))
percent_excl <- as.list((table(d$cond) - table(d_clean$cond))/table(d$cond))

## get mean age and gender:
mean_age = mean(as.numeric(d$age), na.rm = TRUE); mean_age
gender_f = table(d$gender)[2]/sum(table(d$gender)); gender_f

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_merged <- array(dim=c(0, 13))
colnames(d_merged) <- c('vA_sue', 'vB_sue', 'defec', 'negl', 'countf', 'avoid', 'avoid_fault', 
                        'superh', 'comp1', 'comp2', 'mod', 'age', 'cond_name')
d_merged <- as.data.frame(d_merged, stringsAsFactors=FALSE) 

## assess moderator data 
moderator_mat = d_clean[54:58]
moderator_mat$av_trust_5_1 = 100 - as.numeric(moderator_mat$av_trust_5_1) # reverse code one moderator
moderator_mat <- data.frame(sapply(moderator_mat, as.numeric))
cb_alpha = cronbach.alpha(moderator_mat) # calculate cronbach alpha
d_clean$moderator <- rowMeans(moderator_mat) # averaged moderator measure
d_clean <- d_clean %>% relocate(moderator, .after=comp_accident) # move moderator measure

## extract only columns containing measures and other relevant info for analysis
fixed_cols = c(51:54,64,73) # fixed columns - superhuman+comp checks+mod, age, conditions
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$av), ] <- subset(d_clean, cond == "av")[c(29:35,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$hdv), ] <- subset(d_clean, cond == "hdv")[c(44:50,fixed_cols)]

# make columns numeric
d_merged[,1:12] <- lapply(d_merged[,1:12], as.numeric)

## assign trust levels where low trust=1, high trust=2
d_merged$trust_level <- ifelse(d_merged$mod>50, "high", "low")
d_merged$trust_level_n <- ifelse(d_merged$trust_level=="high",2,1)

# cond_n where av=1, human=2
d_merged$cond_n <- ifelse(d_merged$cond_name=="av", 1, 2)

## ================================================================================================================
##                                             DATA ANALYSIS - T-TESTS               
## ================================================================================================================

cor(d_merged[,c(1:8,11)]) # check correlations between measures

## Sue, at-fault (DV)
t.test(vA_sue ~ cond_name, data = d_merged)
aggregate(vA_sue ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$vA_sue, d_merged$cond_name)

## Sue, not-at-fault (DV)
t.test(vB_sue ~ cond_name, data = d_merged)
aggregate(vB_sue ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$vB_sue, d_merged$cond_name)

## Defective (M)
t.test(defec ~ cond_name, data = d_merged)
aggregate(defec ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$defec, d_merged$cond_name)

## Counterfactual (M)
t.test(countf ~ cond_name, data = d_merged)
aggregate(countf ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$countf, d_merged$cond_name)

## Could avoid (M)
t.test(avoid ~ cond_name, data = d_merged)
aggregate(avoid ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$avoid, d_merged$cond_name)

## Trust (MOD)
t.test(mod ~ cond_name, data = d_merged)
aggregate(mod ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$mod, d_merged$cond_name)

## AV superhuman expectations (MOD)
t.test(superh ~ cond_name, data = d_merged)
aggregate(superh ~ cond_name, data = d_merged, FUN = sd)
cohen.d(d_merged$superh, d_merged$cond_name)

## trust agreement with counterfactual
t.test(countf ~ trust_level_n, data = d_merged[d_merged$cond_name=="av", ])
t.test(countf ~ trust_level_n, data = d_merged[d_merged$cond_name=="hdv", ])

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

mediation <- FALSE 

if(mediation) {
  
  source("../../process.R")
  
  # test age as moderator
  summary(lm(vB_sue ~ cond_n*age, data=d_merged))
  
  # SIMPLE MEDIATION
  process(data = d_merged, y = "vB_sue", x = "cond_n", 
          m =c("countf"), model = 4, effsize =1, total =1, stand =1, 
          contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
  
  # MODERATED MEDIATION
  # the effect of trust on A path (7)
  process(data = d_merged, y = "vB_sue", x = "cond_n", 
          m =c("countf"), w = "mod", model = 7, effsize =1, total =1, stand =1, 
          contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
  
  # PARALLEL MEDIATION
  # investigate alternative mediations: (i) trust; (ii) av expectations
  process(data = d_merged, y = "vB_sue", x = "cond_n",
          m =c("mod"), model = 4, effsize =1, total =1, stand =1,
          contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
  
  process(data = d_merged, y = "vB_sue", x = "cond_n",
          m =c("superh"), model = 4, effsize =1, total =1, stand =1,
          contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
  
}

## ================================================================================================================
##                                            PLOTTING COUNTF AGREEMENT                
## ================================================================================================================

t_labels <- c("AV", "HDV")
fill_labels <- c("High trust in AVs", "Low trust in AVs")

## trust v. counterfactual relationship
p_val_L <- t.test(countf ~ trust_level_n, data = d_merged[d_merged$cond_name=="av", ])$p.value
p_val_R <- t.test(countf ~ trust_level_n, data = d_merged[d_merged$cond_name=="hdv", ])$p.value
p1_0 <- plot_2x2(d_merged, x=cond_name, y=countf, fill=trust_level, p_val_L, p_val_R, 
                 title="Agreement Wt. Counterfactual", t_labels, fill_labels)
p1_0 <- p1_0 + xlab ("Vehicle Type") + ylab ("Mean Agreement")
p1_0

png(file = "../plots/s1_trust.png", width = 2*900, height = 2*500, res = 200)  # width and height are in inches
plot(p1_0)
dev.off()

## ================================================================================================================
##                                            PLOTTING MAIN MEASURES                
## ================================================================================================================

t_labels <- c("AV", "HDV")
sig_comparisons <- c("av", "hdv")

## Sue, at-fault (DV)
p_val = t.test(vA_sue ~ cond_name, data = d_merged)$p.value
p2_1 <- plot_std(d_merged, x=cond_name, y=vA_sue, p_val, 
                 title="Sue Veh. A Driver", t_labels, sig_comparisons)

## Sue, not-at-fault (DV)
p_val = t.test(vB_sue ~ cond_name, data = d_merged)$p.value
p2_2 <- plot_std(d_merged, x=cond_name, y=vB_sue, p_val, 
                 title="Sue Veh. B Manufacturer", t_labels, sig_comparisons)

## Counterfactual (M)
p_val = t.test(countf ~ cond_name, data = d_merged)$p.value
p2_3 <- plot_std(d_merged, x=cond_name, y=countf, p_val, 
                 title="Counterfactual", t_labels, sig_comparisons)

## Defective (M)
p_val = t.test(defec ~ cond_name, data = d_merged)$p.value
p2_4 <- plot_std(d_merged, x=cond_name, y=defec, p_val, 
                 title="Veh. B Defective", t_labels, sig_comparisons)

## Could have done more to avoid (M)
p_val = t.test(avoid ~ cond_name, data = d_merged)$p.value
p2_5 <- plot_std(d_merged, x=cond_name, y=avoid, p_val, 
                 title="Could Have Done More", t_labels, sig_comparisons)

figure2 <- ggarrange(p2_1, p2_2, p2_3, p2_4, p2_5, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure2 <- annotate_figure(figure2,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 
figure2

png(file = "../plots/s1_plot.png", width = 2*900, height = 2*700, res = 200)  # width and height are in inches
plot(figure2)
dev.off()

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================
