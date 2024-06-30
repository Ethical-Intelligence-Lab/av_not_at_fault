## clear workspace
rm(list = ls()) 

options(download.file.method="libcurl") # sets method for downloading files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to current directory

source("../common.R") # install packages; import common plotting functions

## ================================================================================================================
##                                            PRE-PROCESSING (EXCLUSIONS)                 
## ================================================================================================================

## read in data: 
d <- read.csv('e3_countf_constr.csv')

## explore dataframe: 
dim(d) # provide dimensions of the dataframe by row [1] and column [2]
colnames(d)
summary(d)

## perform attention exclusions:
# remove responses from the dataframe that failed attention checks (i.e., "1" or "2")
d <- subset(d, (d$att1 == 2 & d$att2 == 2))
dim(d) # number of participants should decrease after attention exclusions

## new columns to label agent (av/hdv), scenario (cnstr/uncnstr), and condition (combination of the 2)
d$agent <- ifelse(d$FL_4_DO %in% c("FL_39", "FL_54"), "av", "hdv") #FL_40/58 are hdv
d$scen <- ifelse(d$FL_4_DO %in% c("FL_39", "FL_40"), "cnstr", "uncnstr") #FL_54/58 are uncnstr
d$cond <- paste(d$agent, d$scen, sep="_")

## get number of participants before comp exclusions
n_orig_all <- dim(d)[1]
n_orig <- as.list(table(d$cond))

## perform comp exclusions
d_clean <- d
d_clean <- subset(d_clean, comp_accident == 1) # include only comp check 2 passes
d_clean <- subset(d_clean, !(agent == "av" & comp1 != 1)) # remove comp check 1 fails for av
d_clean <- subset(d_clean, !(agent == "hdv" & comp1 != 2)) # remove comp check 1 fails for hdv

## get number of participants AFTER conclusions
n_final_all <- dim(d_clean)[1]
percent_excl_all <- (n_orig_all - n_final_all)/n_orig_all
n_excl_all <- n_orig_all - n_final_all
n_final <- as.list(table(d_clean$cond))
percent_excl <- as.list((table(d$cond) - table(d_clean$cond))/table(d$cond))

## duplicate AV condition vB liability columns to match with HDV condition driver liability cols
d_clean$vB_frm_liab_AV_cnstr_2 <- d_clean$vB_frm_liab_AV_cnstr_1 # duplicate
d_clean <- d_clean %>% relocate(vB_frm_liab_AV_cnstr_2, .after=vB_frm_liab_AV_cnstr_1) # move new column
d_clean$vB_frm_liab_AV_2 <- d_clean$vB_frm_liab_AV_1 # duplicate
d_clean <- d_clean %>% relocate(vB_frm_liab_AV_2, .after=vB_frm_liab_AV_1) # move new column

## get mean age and gender:
mean_age = mean(as.numeric(d$age), na.rm = TRUE) # removing NAs from the dataframe before computing mean 
gender_f = table(d$gender)["2"]/sum(table(d$gender))

## ================================================================================================================
##                                                   SUBSETTING                 
## ================================================================================================================

## create df to fill
d_merged <- array(dim=c(0, 11))
colnames(d_merged) <- c('vA_liable', 'vB_m_v_d_liable', 'vB_m_v_m_liable', 'vA_cntrfctl', 'vB_cntrfctl', 
                        'avoid', 'comp1', 'comp2', 'age', 'agent_name', 'scen_name')
d_merged <- as.data.frame(d_merged, stringsAsFactors=FALSE)

## select only the used columns
fixed_cols = c(77:78,82,97:98) # fixed columns - comp checks, age, conditions
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$av_cnstr), ] <- subset(d_clean, cond == "av_cnstr")[c(29:34,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$hdv_cnstr), ] <- subset(d_clean, cond == "hdv_cnstr")[c(43:48,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$av_uncnstr), ] <- subset(d_clean, cond == "av_uncnstr")[c(57:62,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$hdv_uncnstr), ] <- subset(d_clean, cond == "hdv_uncnstr")[c(71:76,fixed_cols)]

# make columns numeric
d_merged[,1:9] <- lapply(d_merged[,1:9], as.numeric)

# agent_n where av=1, human=2; scen_n where cnstr=1, uncnstr=2
d_merged$agent_n <- ifelse(d_merged$agent_name=="av", 1, 2)
d_merged$scen_n <- ifelse(d_merged$scen_name=="cnstr", 1, 2)

## ================================================================================================================
##                                              DATA ANALYSIS - ANOVA              
## ================================================================================================================

cor(d_merged[,1:6]) # check correlations between measures

## Liable, Manufacturer vs Manufacturer (DV)
m_v_m_mod <- aov(vB_m_v_m_liable ~ as.factor(agent_n) * as.factor(scen_n), data = d_merged)
summary(m_v_m_mod)
anova_stats(m_v_m_mod)

# summary statistics
aggregate(vB_m_v_m_liable ~ agent_name, data = d_merged, FUN = mean)
aggregate(vB_m_v_m_liable ~ scen_name, data = d_merged, FUN = mean)

## Liable, Manufacturer vs Driver (DV)
m_v_d_mod <- aov(vB_m_v_d_liable ~ as.factor(agent_n) * as.factor(scen_n), data = d_merged)
summary(m_v_d_mod)
anova_stats(m_v_d_mod)

# summary statistics
aggregate(vB_m_v_d_liable ~ agent_name, data = d_merged, FUN = mean)
aggregate(vB_m_v_d_liable ~ scen_name, data = d_merged, FUN = mean)

## ================================================================================================================
##                                            DATA ANALYSIS - T-TESTS               
## ================================================================================================================

## Liable, at-fault (DV)
t.test(vA_liable ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ])
t.test(vA_liable ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ])
cohen.d(d_merged[d_merged$scen_name=="cnstr", ]$vA_liable, d_merged[d_merged$scen_name=="cnstr", ]$agent_name)
cohen.d(d_merged[d_merged$scen_name=="uncnstr", ]$vA_liable, d_merged[d_merged$scen_name=="uncnstr", ]$agent_name)

## Liable, Manufacturer vs Manufacturer (DV)
t.test(vB_m_v_m_liable ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ])
t.test(vB_m_v_m_liable ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ])
cohen.d(d_merged[d_merged$scen_name=="cnstr", ]$vB_m_v_m_liable, d_merged[d_merged$scen_name=="cnstr", ]$agent_name)
cohen.d(d_merged[d_merged$scen_name=="uncnstr", ]$vB_m_v_m_liable, d_merged[d_merged$scen_name=="uncnstr", ]$agent_name)

## Liable, Manufacturer vs Driver (DV)
t.test(vB_m_v_d_liable ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ])
t.test(vB_m_v_d_liable ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ])
cohen.d(d_merged[d_merged$scen_name=="cnstr", ]$vB_m_v_d_liable, d_merged[d_merged$scen_name=="cnstr", ]$agent_name)
cohen.d(d_merged[d_merged$scen_name=="uncnstr", ]$vB_m_v_d_liable, d_merged[d_merged$scen_name=="uncnstr", ]$agent_name)

## Vehicle A Counterfactual (M)
t.test(vA_cntrfctl ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ])
t.test(vA_cntrfctl ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ])
cohen.d(d_merged[d_merged$scen_name=="cnstr", ]$vA_cntrfctl, d_merged[d_merged$scen_name=="cnstr", ]$agent_name)
cohen.d(d_merged[d_merged$scen_name=="uncnstr", ]$vA_cntrfctl, d_merged[d_merged$scen_name=="uncnstr", ]$agent_name)

## Vehicle B Counterfactual (M)
t.test(vB_cntrfctl ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ])
t.test(vB_cntrfctl ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ])
cohen.d(d_merged[d_merged$scen_name=="cnstr", ]$vB_cntrfctl, d_merged[d_merged$scen_name=="cnstr", ]$agent_name)
cohen.d(d_merged[d_merged$scen_name=="uncnstr", ]$vB_cntrfctl, d_merged[d_merged$scen_name=="uncnstr", ]$agent_name)

# Could have done more to avoid (M)
t.test(avoid ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ])
t.test(avoid ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ])
cohen.d(d_merged[d_merged$scen_name=="cnstr", ]$avoid, d_merged[d_merged$scen_name=="cnstr", ]$agent_name)
cohen.d(d_merged[d_merged$scen_name=="uncnstr", ]$avoid, d_merged[d_merged$scen_name=="uncnstr", ]$agent_name)

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

mediation <- FALSE #change to true if you want to run this code

if(mediation) {
    source("../process.R")
  
    # test age as moderator
    summary(lm(vB_m_v_m_liable ~ agent_n*age, data=d_merged))
    
    # SERIAL MEDIATION
    process(data = d_merged, y = "vB_m_v_m_liable", x = "agent_n", 
            m =c("vB_cntrfctl", "avoid"), model = 6, effsize =1, total =1, stand =1, 
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
    
    # MODERATED SERIAL MEDIATION
    # the effect of scenario on center path (91)
    process(data = d_merged, y = "vB_m_v_m_liable", x = "agent_n", 
            m =c("vB_cntrfctl", "avoid"), w = "scen_n", model = 91, effsize =1, total =1, stand =1, 
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
}

## ================================================================================================================
##                                              PLOTTING 2X2 FIGURE                 
## ================================================================================================================

t_labels <- c("Constrained", "Unconstrained")
fill_labels <- c("AV", "HDV")

## Liable, Manufacturer vs Manufacturer (DV)
p_val_L = t.test(vB_m_v_m_liable ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ], paired = FALSE)$p.value
p_val_R = t.test(vB_m_v_m_liable ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ], paired = FALSE)$p.value
p1_1 <- plot_2x2(d_merged, x=scen_name, y=vB_m_v_m_liable, fill=agent_name, p_val_L, p_val_R, 
                 title="Veh. B Manufacturer Liability", t_labels, fill_labels)

## Liable, Manufacturer vs Driver (DV)
p_val_L = t.test(vB_m_v_d_liable ~ agent_name, data = d_merged[d_merged$scen_name=="cnstr", ], paired = FALSE)$p.value
p_val_R = t.test(vB_m_v_d_liable ~ agent_name, data = d_merged[d_merged$scen_name=="uncnstr", ], paired = FALSE)$p.value
p1_2 <- plot_2x2(d_merged, x=scen_name, y=vB_m_v_d_liable, fill=agent_name, p_val_L, p_val_R, 
                 title="Veh. B Manufacturer\nor Driver Liability", t_labels, fill_labels)

figure1 <- ggarrange(p1_1, p1_2, nrow=1, ncol=2,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure1 <- annotate_figure(figure1, left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Scenario Type", color="black", face ="plain",size=18)) 

png(file = "../plots/e3_constraints.png", width = 2*900, height = 2*500, res = 200)  # width and height are in inches
plot(figure1)
dev.off()

## ================================================================================================================
##                                              PLOTTING BY AGENT               
## ================================================================================================================

t_labels <- c("AV", "HDV")
sig_comparisons <- c("av", "hdv")

## Liable, at-fault (DV)
p_val = t.test(vA_liable ~ agent_name, data = d_merged, paired = FALSE)$p.value
p2_1 <- plot_std(d_merged, x=agent_name, y=vA_liable, p_val, 
                 title="Veh. A Driver Liability", t_labels, sig_comparisons)

## Liable, Manufacturer vs Manufacturer (DV)
p_val = t.test(vB_m_v_m_liable ~ agent_name, data = d_merged, paired = FALSE)$p.value
p2_2 <- plot_std(d_merged, x=agent_name, y=vB_m_v_m_liable, p_val, 
             title="Veh. B Manufacturer Liability", t_labels, sig_comparisons)

## Liable, Manufacturer vs Driver (DV)
p_val = t.test(vB_m_v_d_liable ~ agent_name, data = d_merged, paired = FALSE)$p.value
p2_3 <- plot_std(d_merged, x=agent_name, y=vB_m_v_d_liable, p_val, 
             title="Veh. B Manufacturer\nor Human Driver Liability", t_labels, sig_comparisons)

## Vehicle A Counterfactual (M)
p_val = t.test(vA_cntrfctl ~ agent_name, data = d_merged, paired = FALSE)$p.value
p2_4 <- plot_std(d_merged, x=agent_name, y=vA_cntrfctl, p_val, 
             title="Consider Veh. A Counterfactual", t_labels, sig_comparisons)

## Vehicle B Counterfactual (M)
p_val = t.test(vB_cntrfctl ~ agent_name, data = d_merged, paired = FALSE)$p.value
p2_5 <- plot_std(d_merged, x=agent_name, y=vB_cntrfctl, p_val, 
             title="Consider Veh. B Counterfactual", t_labels, sig_comparisons)

## Could have done more to avoid (M)
p_val = t.test(avoid ~ agent_name, data = d_merged, paired = FALSE)$p.value
p2_6 <- plot_std(d_merged, x=agent_name, y=avoid, p_val, 
             title="Could Have Done More", t_labels, sig_comparisons)

figure2 <- ggarrange(p2_1, p2_2, p2_3, p2_4, p2_5, p2_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure2 <- annotate_figure(figure2,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 

png(file = "../plots/e3_plot.png", width = 2*900, height = 2*700, res = 200)  # width and height are in inches
plot(figure2)
dev.off()

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================