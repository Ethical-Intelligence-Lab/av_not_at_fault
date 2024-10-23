### common install / functions required for av_not_at_fault

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
               'Hmisc',           # confidence intervals
               'sjstats',         # anova stats
               'lavaan',          # cfa tests for discriminant validity
               'semTools'         # discriminant validity tests       
               
)

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

# plotting function
plot_std <- function(d_merged, x, y, p_val, title, t_labels, sig_comparisons) {

  anno <- get_annotation(p_val)
  
  p <- ggplot(d_merged,aes(x=factor({{x}}),y={{y}})) +  
    theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+ 
    geom_signif(comparisons = list(sig_comparisons), annotation=unlist(anno[1]), textsize = unlist(anno[2]))
  
  p <- p + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_x_discrete(labels=t_labels) +
    ggtitle(title) +
    xlab ("") + ylab ("") +
    theme_classic() +
    theme(axis.text.x = element_text(size=12)) +
    theme(axis.text.y = element_text(size=10)) +
    theme(plot.title = element_text(size=12, hjust=0.5)) +
    geom_bar(stat="summary", position = position_dodge(), width = 0.9, alpha = 0.38, size = 0.75) +
    # geom_violin(width=0.9, alpha=0.38, size=0.75) +  
    # geom_sina(alpha=0.6, size=0.95, color = "#999999") +
    stat_summary(fun.data = "mean_cl_boot", color = "black", 
                 size=0.4, 
                 position = position_dodge(width = 0.9)) +
    stat_summary(fun.data = "mean_cl_boot", color = "black", 
                 position = position_dodge(width = 0.9),
                 geom="errorbar", width = 0.2)
  return(p)

}

## plotting function (2x2)
plot_2x2 <- function(d_merged, x, y, fill, p_val_L, p_val_R, title, t_labels, fill_labels) {
  
  L_anno <- get_annotation(p_val_L)
  R_anno <- get_annotation(p_val_R)
  
  p <- ggplot(d_merged, aes(x=factor({{x}}), y={{y}}, fill={{fill}})) +  
    theme_bw() + coord_cartesian(ylim=c(1,110))+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
    geom_signif(y_position = 105.00, xmin = c(0.8,1.8), xmax = c(1.2,2.2), annotation = c(unlist(L_anno[1]),unlist(R_anno[1])), textsize=7.5)
  
  p <- p + theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_x_discrete(labels=t_labels) +
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
    geom_bar(stat="summary", position = position_dodge(), width = 0.9, alpha = 0.38, size = 0.75) +
    stat_summary(fun.data = "mean_cl_boot", color = "black", 
                 size=0.4, 
                 position = position_dodge(width = 0.9)) +
    stat_summary(fun.data = "mean_cl_boot", color = "black", 
                 position = position_dodge(width = 0.9),
                 geom="errorbar", width = 0.2)
  return(p)
  
}