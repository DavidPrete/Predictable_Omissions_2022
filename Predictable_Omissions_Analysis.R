#library(plyr)
#library(ggpubr)
#library(emmeans)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(colorspace)

rm(list = ls())
#DIR from home computer
#setwd("C:/Users/its_d/Dropbox/Predictability and Omission/Analysis/")

# #DIR from office computer
setwd("D:/David/DropBox/Dropbox/Predictability and Omission/Analysis")

dir()

#----- Helper Functions ----

box_plot_amps = function(data, ERP_comp, title){
  
  if(missing(title)){
    title = ERP_comp
  }
  
  save_name = paste0(ERP_comp,'.tiff')
  
  hex_colors = c("#66B6E2","#B11818")
    
  bxplt = ggplot(data, aes(x = Predictability, y = Mean_Amp, fill = Predictability)) +
    geom_boxplot(color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.15,size =0.75)+
    labs( y = expression(paste("Amplitude ( ", mu, "v)")))+
    facet_grid(cols = vars(Laterality), rows = vars(Centrality))+
    ggtitle(title)+
    scale_fill_manual(values = hex_colors)+
    # scale_color_manual(values = darken(hex_colors, amount = 0.5))+
    theme_bw()+
    theme(panel.background = element_blank(),
            axis.text.y  = element_text(size=13),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=17),
            strip.text.x = element_text(size=15),
            strip.text.y = element_text(size=15),
            plot.title   = element_text(size=20, hjust = 0.5), 
            legend.title = element_blank(),
            legend.text  = element_text(size=15),
            legend.key.size = unit(1.5, 'cm'),
            legend.position="bottom"
            )
  bxplt

  ggsave(save_name,units="mm", width=115, height=90, dpi=600)
  bxplt
  
}

rep_anovas = function(data, effect_size){
  
  res.aov <- anova_test(
    data = data, dv = Mean_Amp, wid = Participants,
    within = c(Predictability,Centrality, Laterality),
    effect.size =effect_size, 
    type =2
  )
  get_anova_table(res.aov)
}

post_hoc_tests = function(data, check){
  
  if(check==0){
    dfs = split(data, f = list(data$Predictability,data$Laterality))
    conditions = names(dfs)
  }else{
    data = data %>% 
      group_by(Participants,Predictability,Laterality) %>%
      get_summary_stats(Mean_Amp, type = "mean_sd") 
    data = rename(data, Mean_Amp = mean)
    dfs = split(data, f = list(data$Predictability,data$Laterality))
    
  }
  
  count = 1 
  ttest_lists = vector(mode='list', length=length(dfs))
  diffs_list  = vector(mode = 'list',length = length(dfs)/2)

  for(ii in seq(from=1, to=length(dfs), by=2)){
    jj = ii +1
    pred = dfs[[conditions[ii]]]$Mean_Amp
    unpred = dfs[[conditions[jj]]]$Mean_Amp
    diff = unpred - pred
    
    t_tests = t.test(unpred,pred , paired = TRUE)
    ttest_lists[[count]] = t_tests
    diffs_list[[count]]  = diff
    count = count+ 1
  }
  
  left_diffs   = diffs_list[[1]] 
  middle_diffs = diffs_list[[2]]
  right_diffs  = diffs_list[[3]] 
  
  #testing diff scores from left region and middle regions
  ttest_lists[[4]] = t.test(left_diffs,middle_diffs,paired = TRUE)
  #testing diff scores from the right region and the middles regions
  ttest_lists[[5]] = t.test(right_diffs,middle_diffs,paired = TRUE)
  #testing the diff scores from the left region and right region
  ttest_lists[[6]] = t.test(left_diffs,right_diffs,paired = TRUE)
  
    
  
  return(ttest_lists)
}

test_against_zero = function(data){
  ttest_list<- vector(mode='list', length=2)
  
  data = data %>% 
    group_by(Participants,Predictability) %>%
    get_summary_stats(Mean_Amp, type = "mean_sd") 
  
  
  dfs = split(data, f = list(data$Predictability))
  pred_comp   = dfs[[1]] 
  unpred_comp = dfs[[2]]
  
  pred_test   = t.test(pred_comp$mean,mu = 0)
  unpred_test = t.test(unpred_comp$mean,mu = 0)
  ttest_list[[1]] = pred_test
  ttest_list[[2]] = unpred_test
  
  return(ttest_list)
  
}

post_hoc_aov = function(data, check){
  
  front = data %>% filter(Centrality =="Frontal") %>% subset(select = -c(Centrality))
  centr = data %>% filter(Centrality =="Central") %>% subset(select = -c(Centrality ))
  
  front_aov <- anova_test(
    data = front, dv = Mean_Amp, wid = Participants,
    within = c(Predictability,Laterality),
    effect.size =effect_size
    )
  
  centr_aov <- anova_test(
    data = centr, dv = Mean_Amp, wid = Participants,
    within = c(Predictability,Laterality),
    effect.size =effect_size
    )

  if(front_aov$ANOVA$p[3]<0.05){
    front_posts = post_hoc_tests(front,0)
  }else{
    front_posts = list()
  }
  
  if(centr_aov$ANOVA$p[3]<0.05){
    centr_posts = post_hoc_tests(centr,0)
  }else{
    centr_posts = list()
  }
 
  anovas = list(front_aov, centr_aov,front_posts,centr_posts)
  return(anovas)
}


#----- Load and Convert data -----

# READ IN MEAN AMP DATA 
data = read.csv("Predictable_Omissions_Mean_Amps2.csv")
colnames(data)

data_long = data %>% pivot_longer(
  cols = names(data[-1]),
  names_to = c("Component","Predictability","Centrality","Laterality"),
  names_sep="_",
  values_to="Mean_Amp"
)

data_long$Centrality     = factor(data_long$Centrality, levels = c("Frontal","Central"))
data_long$Laterality     = factor(data_long$Laterality, levels = c("Left", "Middle", "Right"))
data_long$Predictability = factor(data_long$Predictability, levels = c("Predictable", "Unpredictable"))


P1_data = data_long %>% filter(Component == "P1") %>% subset(select = -c(Component))
# N1_data = data_long %>% filter(Component == "N1" | Component == "Baseline")
# P2_data = data_long %>% filter(Component == "P2" | Component == "Baseline")
# N2_data = data_long %>% filter(Component == "N2" | Component == "Baseline")

N1Diff_data =  data_long %>% filter(Component == "N1Diff") %>% subset(select = -c(Component))
P2Diff_data =  data_long %>% filter(Component == "P2Diff") %>% subset(select = -c(Component))
N2Diff_data =  data_long %>% filter(Component == "N2Diff") %>% subset(select = -c(Component))

# N1Diff_data2 =  data_long %>% filter(Component == "N1Diff"| Component == "Baseline")
# P2Diff_data2 =  data_long %>% filter(Component == "P2Diff"| Component == "Baseline")
# N2Diff_data2 =  data_long %>% filter(Component == "N2Diff"| Component == "Baseline")

#----- Visualize Data ----

#setwd("C:/Users/David/Dropbox/Predictability and Omission/figures/")
setwd("D:/David/DropBox/Dropbox/Predictability and Omission/figures/")

box_plot_amps(P1_data, "P1 Mean Amplitude")
# box_plot_amps(P2_data, "P2 Mean Amplitude")
# box_plot_amps(N2_data, "N2 Mean Amplitude")

box_plot_amps(N1Diff_data,"N1-P1 Mean Amplitude",expression(paste(N1[d]," Mean Amplitude")))
box_plot_amps(P2Diff_data,"P2-N1 Mean Amplitude",expression(paste(P2[d]," Mean Amplitude")))
box_plot_amps(N2Diff_data,"N2-P2 Mean Amplitude",expression(paste(N2[d]," Mean Amplitude")))


#----- ANOVAs ----
#P1 ANVOA
effect_size = "pes"

#Repeated Measures ANOVA for the individual mean amplitude for each ERP componenet
rep_anovas(P1_data, effect_size)
rep_anovas(N1_data, effect_size)
rep_anovas(P2_data, effect_size)
rep_anovas(N2_data, effect_size)

#Repeated Measures ANOVA for the individual mean amplitude for each ERP componenet relative to the previous component
rep_anovas(N1Diff_data, effect_size)
rep_anovas(P2Diff_data, effect_size)
rep_anovas(N2Diff_data, effect_size)

#----- Post-Hoc Tests -----
 
P1_post_hocs = post_hoc_aov(P1_data)
p1_vs_base = test_against_zero(P1_data)
# N1_post_hocs     = post_hoc_aov(N1_data)
# P2_post_hocs     = post_hoc_aov(P2_data)

N1Diff_post_hocs = post_hoc_aov(N1Diff_data)
P2Diff_post_hocs = post_hoc_aov(P2Diff_data)
#N2Diff_post_hocs = post_hoc_aov(N2Diff_data)
N2Diff_post_hocs = post_hoc_tests(N2Diff_data,1)


# testing to see if the ERP components are significantly different from baseline 
n1_vs_base = test_against_baseline(N1Diff_data2)
p2_vs_base = test_against_baseline(P2Diff_data2)
n2_vs_base = test_against_baseline(N2Diff_data2)

# testing to see if the ERP components are significantly different from hypothetical zero 
n1_vs_base = test_against_zero(N1Diff_data)
p2_vs_base = test_against_zero(P2Diff_data)
n2_vs_base = test_against_zero(N2Diff_data)





