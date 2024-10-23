## clear workspace
rm(list = ls()) 

options(download.file.method="libcurl") # sets method for downloading files

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to current directory

source("../common.R") # install packages; import common plotting functions

install.packages("reshape2")
library(reshape2)

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
d <- read.csv('LLM_measures_300_4o_constr.csv')

## explore dataframe: 
dim(d) # will provide dimensions of the dataframe by row [1] and column [2]
colnames(d)
summary(d)

cor(d[,3:4]) # check correlations between measures

d_long <- reshape2::melt(d[, 3:4])

ggplot(d_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplots for two columns", x = "Columns", y = "Values")

## ================================================================================================================
##                                                  PRE-PROCESSING                 
## ================================================================================================================

## read in data: 
d <- read.csv('LLM_measures_300_4o_julia.csv')

## explore dataframe: 
dim(d) # will provide dimensions of the dataframe by row [1] and column [2]
colnames(d)
summary(d)

d_0 <- d[d$scenario == 0, ]
d_1 <- d[d$scenario == 1, ]

cor(d_0[,4:5]) # check correlations between measures
cor(d_1[,4:5]) # check correlations between measures

d0_long <- reshape2::melt(d_0[, 3:4])
d1_long <- reshape2::melt(d_1[, 3:4])

ggplot(d0_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplots for two columns", x = "Columns", y = "Values")

ggplot(d1_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplots for two columns", x = "Columns", y = "Values")

## ================================================================================================================
##                                                    SUBSETTING                 
## ================================================================================================================

## define new data frame to extract pre-processed data into:
d_merged <- array(dim=c(0, 11))
colnames(d_merged) <- c('vA_sue', 'vB_m_v_d_sue', 'vB_m_v_m_sue', 'vA_cntrfctl', 'vB_cntrfctl', 
                        'avoid', 'comp1', 'comp2', 'age', 'agent_name', 'intv_appld')
d_merged <- as.data.frame(d_merged, stringsAsFactors=FALSE) 

## select only the used columns
fixed_cols = c(79:80,84,102:103) # fixed columns - comp checks, age, conditions
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$av_yes), ] <- subset(d_clean, cond == "av_yes")[c(44:49,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$hdv_yes), ] <- subset(d_clean, cond == "hdv_yes")[c(73:78,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$av_no), ] <- subset(d_clean, cond == "av_no")[c(29:34,fixed_cols)]
d_merged[(dim(d_merged)[1]+1):(dim(d_merged)[1]+n_final$hdv_no), ] <- subset(d_clean, cond == "hdv_no")[c(58:63,fixed_cols)]

# make columns numeric
d_merged[,1:9] <- lapply(d_merged[,1:9], as.numeric)

# average Veh. B counterfactual, and Veh. B could have avoided
d_merged$med <- (d_merged$vB_cntrfctl + d_merged$avoid) / 2
d_merged <- d_merged %>% relocate(med, .after=avoid)

# agent_n where av=1, human=2; intv_n where yes=1, no=2
d_merged$agent_n <- ifelse(d_merged$agent_name=="av", 1, 2)
d_merged$intv_n <- ifelse(d_merged$intv_appld=="yes", 1, 2)

## ================================================================================================================
##                                         DATA ANALYSIS - ANOVA              
## ================================================================================================================

cor(d_merged[,1:7]) # check correlations between measures

## Check for discriminant validity
fit.model <- ' M1 =~ vB_cntrfctl
              M2 =~ avoid '
fit <- cfa(fit.model, data = d_merged)
discriminantValidity(fit)

reliability(fit)

## Sue, Manufacturer vs Manufacturer (DV)
m_v_m_mod <- aov(vB_m_v_m_sue ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(m_v_m_mod)
anova_stats(m_v_m_mod)

# summary statistics
aggregate(vB_m_v_m_sue ~ agent_name, data = d_merged, FUN = mean)
aggregate(vB_m_v_m_sue ~ intv_appld, data = d_merged, FUN = mean)

## Sue, Manufacturer vs Driver (DV)
m_v_d_mod <- aov(vB_m_v_d_sue ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(m_v_d_mod)
anova_stats(m_v_d_mod)

# summary statistics
aggregate(vB_m_v_d_sue ~ agent_name, data = d_merged, FUN = mean)
aggregate(vB_m_v_d_sue ~ intv_appld, data = d_merged, FUN = mean)

## Vehicle B Counterfactual (M)
vB_countf_mod <- aov(vB_cntrfctl ~ as.factor(agent_n) * as.factor(intv_n), data = d_merged)
summary(vB_countf_mod)
anova_stats(vB_countf_mod)

# summary statistics
aggregate(vB_cntrfctl ~ agent_name, data = d_merged, FUN = mean)
aggregate(vB_cntrfctl ~ intv_appld, data = d_merged, FUN = mean)

## ================================================================================================================
##                                            DATA ANALYSIS - T-TESTS               
## ================================================================================================================

## Sue, at-fault (DV)
t.test(vA_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(vA_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$vA_sue, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$vA_sue, d_merged[d_merged$intv_appld=="no", ]$agent_name)

## Sue, Manufacturer vs Manufacturer (DV)
t.test(vB_m_v_m_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(vB_m_v_m_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$vB_m_v_m_sue, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$vB_m_v_m_sue, d_merged[d_merged$intv_appld=="no", ]$agent_name)

## Sue, Manufacturer vs Driver (DV)
t.test(vB_m_v_d_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(vB_m_v_d_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$vB_m_v_d_sue, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$vB_m_v_d_sue, d_merged[d_merged$intv_appld=="no", ]$agent_name)

## Vehicle A Counterfactual (M)
t.test(vA_cntrfctl ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(vA_cntrfctl ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$vA_cntrfctl, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$vA_cntrfctl, d_merged[d_merged$intv_appld=="no", ]$agent_name)

## Vehicle B Counterfactual (M)
t.test(vB_cntrfctl ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(vB_cntrfctl ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$vB_cntrfctl, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$vB_cntrfctl, d_merged[d_merged$intv_appld=="no", ]$agent_name)

# Could have done more to avoid (M)
t.test(avoid ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(avoid ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$avoid, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$avoid, d_merged[d_merged$intv_appld=="no", ]$agent_name)

# Averaged mediator (M)
t.test(med ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])
t.test(med ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])
cohen.d(d_merged[d_merged$intv_appld=="yes", ]$med, d_merged[d_merged$intv_appld=="yes", ]$agent_name)
cohen.d(d_merged[d_merged$intv_appld=="no", ]$med, d_merged[d_merged$intv_appld=="no", ]$agent_name)

## ================================================================================================================
##                                             MEDIATION ANALYSIS              
## ================================================================================================================

mediation <- FALSE #change to true if you want to run this code

if(mediation) {
    source("../process.R")
    
    # test age as moderator
    summary(lm(vB_m_v_m_sue ~ agent_n*age, data=d_merged))
    
    # SERIAL MEDIATION
    process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n",
            m =c("vB_cntrfctl", "avoid"), model = 6, effsize =1, total =1, stand =1,
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
    
    # SERIAL MEDIATION (flipped)
    process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n", 
            m =c("avoid", "vB_cntrfctl"), model = 6, effsize =1, total =1, stand =1, 
            contrast =1, boot = 10000 , modelbt = 1, seed = 654322)
    
    # MODERATED SERIAL MEDIATION
    # the effect of intervention on A path (83)
    process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n",
            m =c("vB_cntrfctl", "avoid"), w = "intv_n", model = 83, effsize =1, total =1, stand =1,
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
    
    # MODERATED SERIAL MEDIATION
    # the effect of intervention on center path (91)
    process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n",
            m =c("vB_cntrfctl", "avoid"), w = "intv_n", model = 91, effsize =1, total =1, stand =1,
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
    
    # PARALLEL MEDIATION (averaged mediators)
    process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n", 
            m =c("med"), model = 4, effsize =1, total =1, stand =1, 
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
    
    # MODERATED MEDIATION (averaged mediators)
    # the effect of intervention on A path (7)
    process(data = d_merged, y = "vB_m_v_m_sue", x = "agent_n", 
            m =c("med"), w = "intv_n", model = 7, effsize =1, total =1, stand =1, 
            contrast =1, boot = 10000 , modelbt = 1, seed = 654321)
}

## ================================================================================================================
##                                              PLOTTING 2X2 FIGURE                 
## ================================================================================================================

t_labels <- c("No Intervention", "Intervention")
fill_labels <- c("AV", "HDV")

## Sue, Manufacturer vs Manufacturer (DV)
p_val_L = t.test(vB_m_v_m_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])$p.value
p_val_R = t.test(vB_m_v_m_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])$p.value
p1_1 <- plot_2x2(d_merged, x=intv_appld, y=vB_m_v_m_sue, fill=agent_name, p_val_L, p_val_R, 
                 title="Sue Veh. B Manufacturer", t_labels, fill_labels)

## Sue, Manufacturer vs Manufacturer (DV)
p_val_L = t.test(vB_m_v_d_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])$p.value
p_val_R = t.test(vB_m_v_d_sue ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])$p.value
p1_2 <- plot_2x2(d_merged, x=intv_appld, y=vB_m_v_d_sue, fill=agent_name, p_val_L, p_val_R, 
                 title="Sue Veh. B Manufacturer\nor Human Driver", t_labels, fill_labels)

## Vehicle B Counterfactual (M)
p_val_L = t.test(vB_cntrfctl ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])$p.value
p_val_R = t.test(vB_cntrfctl ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])$p.value
p1_3 <- plot_2x2(d_merged, x=intv_appld, y=vB_cntrfctl, fill=agent_name, p_val_L, p_val_R, 
                 title="Consider Veh. B Counterfactual", t_labels, fill_labels)

## Could have done more to avoid (M)
p_val_L = t.test(avoid ~ agent_name, data = d_merged[d_merged$intv_appld=="no", ])$p.value
p_val_R = t.test(avoid ~ agent_name, data = d_merged[d_merged$intv_appld=="yes", ])$p.value
p1_4 <- plot_2x2(d_merged, x=intv_appld, y=avoid, fill=agent_name, p_val_L, p_val_R, 
                 title="Could Have Done More", t_labels, fill_labels)

figure1 <- ggarrange(p1_1, p1_2, p1_3, p1_4, nrow=2, ncol=2,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure1 <- annotate_figure(figure1, left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Intervention Applied", color="black", face ="plain",size=18)) 

png(file = "../plots/e4_intervention.png", width = 2*900, height = 2*900, res = 200)  # width and height are in inches
plot(figure1)
dev.off()

## ================================================================================================================
##                                              PLOTTING BY AGENT               
## ================================================================================================================

t_labels <- c("AV", "HDV")
sig_comparisons <- c("av", "hdv")

## Sue, at-fault (DV)
p_val = t.test(vA_sue ~ agent_name, data = d_merged)$p.value
p2_1 <- plot_std(d_merged, x=agent_name, y=vA_sue, p_val, 
                 title="Sue Veh. A Driver", t_labels, sig_comparisons)

## Liable, Manufacturer vs Manufacturer (DV)
p_val = t.test(vB_m_v_m_sue ~ agent_name, data = d_merged)$p.value
p2_2 <- plot_std(d_merged, x=agent_name, y=vB_m_v_m_sue, p_val, 
                 title="Sue Veh. B Manufacturer", t_labels, sig_comparisons)

## Liable, Manufacturer vs Driver (DV)
p_val = t.test(vB_m_v_d_sue ~ agent_name, data = d_merged)$p.value
p2_3 <- plot_std(d_merged, x=agent_name, y=vB_m_v_d_sue, p_val, 
                 title="Sue Veh. B Manufacturer\nor Human Driver", t_labels, sig_comparisons)

## Vehicle A Counterfactual (M)
p_val = t.test(vA_cntrfctl ~ agent_name, data = d_merged)$p.value
p2_4 <- plot_std(d_merged, x=agent_name, y=vA_cntrfctl, p_val, 
                 title="Consider Veh. A Counterfactual", t_labels, sig_comparisons)

## Vehicle B Counterfactual (M)
p_val = t.test(vB_cntrfctl ~ agent_name, data = d_merged)$p.value
p2_5 <- plot_std(d_merged, x=agent_name, y=vB_cntrfctl, p_val, 
                 title="Consider Veh. B Counterfactual", t_labels, sig_comparisons)

## Could have done more to avoid (M)
p_val = t.test(avoid ~ agent_name, data = d_merged)$p.value
p2_6 <- plot_std(d_merged, x=agent_name, y=avoid, p_val, 
                 title="Could Have Done More", t_labels, sig_comparisons)

figure2 <- ggarrange(p2_1, p2_2, p2_3, p2_4, p2_5, p2_6, nrow=2,ncol=3,common.legend = TRUE, legend="top", vjust = 1.0, hjust=0.5) 
figure2 <- annotate_figure(figure2,left = text_grob("Mean Agreement", color="black", face ="plain",size=16, rot=90),
                           bottom = text_grob("Vehicle Type", color="black", face ="plain",size=16)) 

png(file = "../plots/e4_plot.png", width = 2*900, height = 2*700, res = 200)  # width and height are in inches
plot(figure2)
dev.off()

## ================================================================================================================
##                                                  END OF ANALYSIS                 
## ================================================================================================================