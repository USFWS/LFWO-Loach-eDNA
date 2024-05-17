# Purpose:
#   artemis modeling USFWS San Luis Wildlife Refuge Loach livecar 2023

# Code author: 
# Katie Karpenko

# Last updated: 
# March 7, 2023

# Libraries
library(artemis)
library(openxlsx)
library(rdryad)
library(ggplot2)
library(stringr)


# Read data
# dat = read.xlsx("Analysis/Livecar/USFWS_Loach removal_eDNA calibration.xlsx", sheet=1) #update path to data file location
dat = dryad_files_download('xxx')

# set distances plotting order
table(dat$dist_lc_m)
dat$dist_lc_m = factor(dat$dist_lc_m, 
                       level=c("Pre", "10", "50","100","250","500",
                               "Field Control", "Extraction Control"))

# subset data with only targets that were in LC
lc_only=subset(dat, Target !="LSL/FSL")

# Add alpha and beta to df (this value is determined by the standard curve for 
#    each assay). Alpha and beta used in artemis model.
Target = c("BHOO", "ANCH") #target names
alpha = c("22.686", "18.475") #ln(y-intercept) of standard curve for each target name
beta = c("-1.406", "-1.545") #ln(slope) of standard curve for each target name

# create target, alpha, and beta dataframe
ab=data.frame(Target, alpha, beta)

# add alpha and beta to livecar qPCR data
lc_only_alpha_beta=merge(lc_only,ab, by="Target")
lc_only_alpha_beta$alpha=as.numeric(lc_only_alpha_beta$alpha)
lc_only_alpha_beta$beta=as.numeric(lc_only_alpha_beta$beta)

# subset data by canal size
small_canal=subset(lc_only_alpha_beta, location=="Small Channel" & 
                     control!="Control" & dist_lc_m!="Pre")

lrg_canal=subset(lc_only_alpha_beta, location=="Large Channel" & 
                   control!="Control" & dist_lc_m!="Pre")

# Modeling small channel data
m_sc = eDNA_lmer(Cq ~ as.factor(dist_lc_m) + Target + (1|FilterID),
                   data = small_canal,
                   std_curve_alpha = small_canal$alpha,
                   std_curve_beta = small_canal$beta)


m_sc
plot(m_sc)

# export model
saveRDS(m_sc, "Analysis/Livecar/small_channel.rds") # export model results
m_sc = readRDS("Analysis/Livecar/small_channel.rds") # import model results (if already exported)

# Modeling large channel data
m_lc = eDNA_lmer(Cq ~ dist_lc_m + Target + (1|FilterID),
                 data = lrg_canal,
                 std_curve_alpha = lrg_canal$alpha,
                 std_curve_beta = lrg_canal$beta)

# failed. One positive tech rep not enough to model..


# format model output of small canal model
sm2 = as.data.frame(summary(m_sc))
sm2$Predictor = row.names(sm2)
colnames(sm2) = c("mean", "lwr", "median", "upr", "Predictor")
row2 = c("Intercept","50m","100m","250m","500m","Target(Ballyhoo)",
          "sigma ln(eDNA)")
sm2=cbind(sm2,row2)
sm2$row2 = factor(sm2$row2, levels=c("500m", "250m","100m","50m",
                                     "Target(Ballyhoo)","sigma ln(eDNA)", 
                                     "Intercept"))

# plot model output for small canal
sc_model_plot =
  ggplot(sm2[sm2$row2 != "Intercept", ], 
       aes(x = median, 
           y = row2)) +
  geom_point(size = 2.5, color = "gray40") +
  geom_segment(aes(xend = lwr, x = median, y = row2, yend = row2)) +
  geom_segment(aes(xend = upr, x = median, y = row2, yend = row2)) +
  scale_y_discrete(limits=rev)+
  labs(x = "Model Estimate (ln[eDNA])", y = NULL)+
  theme_bw()+
  theme(text=element_text(size=12, family="serif"), 
        plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5))
 

# Export small channel model output
ggsave("Analysis/Manuscript Plots and Figs/sc_model_plot2.jpeg",
       plot=sc_model_plot, 
       width=7, 
       height=1.5*aspect_ratio) # update directory as needed


# p-detect in small channel
# I ran through each combo of factors and manually created table in excel with the output
est_p_detect(variable_levels = c(dist_lc_m50 = 1,
                                 dist_lc_m100= 0,
                                 dist_lc_m250= 0,
                                 dist_lc_m500= 0,
                                 TargetBHOO = 0),
             std_curve_alpha = 18.475, # for loach, from PL_MM assay standard curve
             std_curve_beta = -1.545,  # for loach, from PL_MM assay standard curve
             model_fit = m_sc, n_rep = 3) # 3 reps = 3 filters

# When all combinations are tested with p_detect, read in the manually created PoD excel table
small_canal_pod = read.xlsx("Analysis/Livecar/small_channel_pod.xlsx") # update directory

# plot PoD
sc_pod_plot= 
  ggplot(subset(small_canal_pod, n_techrep=="3"), 
         aes(x = Distance, y = PoD_mean, shape=Target))+
  geom_line()+
  geom_point(size=2)+
  scale_shape_manual(values=c(1,15), name="Estimated Biomass",
                     labels=c("Appx. 16g", "Appx. 100g"),
                     guide = guide_legend(label.theme = element_text(size = 10, 
                                                                     family="serif")))+
  labs(y="Probability of Detection", 
       x="Distance (m) from Fixed Biomass",
       subtitle="Small Canal, Number of filters = 3")+
  theme_bw()+
  theme(text = element_text (size=11, family= "serif"))+
  theme(plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5), legend.position = "bottom")

# export PoD plot
ggsave("Analysis/Manuscript Plots and Figs/sc_pod_plot.jpeg",plot=sc_pod_plot, width=5, 
       height=1.25*aspect_ratio)  # update directory

# PoD plot with error bars
facet_names = as_labeller(
  c(`Anchovy` = "Appx. 16g", `Ballyhoo` = "Appx. 100g"))

sc_pod_plot_eb= 
  ggplot(subset(small_canal_pod, n_techrep=="3"), 
         aes(x = Distance, y = PoD_mean, shape=Target))+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=PoD_low, ymax=PoD_high), width=4) +
  scale_shape_manual(values=c(1,15), 
                     name="Estimated Biomass",
                     labels=c("Appx. 16g", "Appx. 100g"),
                     guide = guide_legend(label.theme = 
                                            element_text(size = 10, 
                                                         family="serif")))+
  labs(y="Probability of Detection", 
       x="Distance (m) from Fixed Biomass")+
  theme_bw()+
  theme(text = element_text (size=12, family= "serif"))+
  theme(plot.title = element_text(hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5), legend.position = "bottom")+
  facet_wrap(~Target, ncol= 2, labeller=facet_names)


ggsave("Analysis/Manuscript Plots and Figs/sc_pod_plot_errorbars.jpeg",
       plot=sc_pod_plot_eb, 
       width=10, 
       height=2*aspect_ratio) # update directory

### END ###
