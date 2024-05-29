# Purpose:
#   artemis modeling USFWS San Luis Wildlife Refuge Loach livecar 2023

# Code author: 
# Katie Karpenko

# Last updated: 
# March 7, 2023

# Libraries
# 
# artemis requires cmdstanr, which you might need to install
#     artemis installation instructions: https://github.com/fishsciences/artemis
#     cmdstanr instructions: https://mc-stan.org/cmdstanr/articles/cmdstanr.html
#         - requires cmdstan and rtools
library(cmdstanr)
# verify path and version
cmdstan_path()
cmdstan_version() # script tested with version "2.34.1"

# If artemis needs to be installed:
# devtools::install_github("fishsciences/artemis") # install artemis
# artemis::compile_models() # don't skip this after successfully installing artemis

library(artemis)
library(openxlsx)
library(rdryad)
library(ggplot2)
library(stringr)


# Read data
# dat = read.xlsx("Analysis/Livecar/USFWS_Loach removal_eDNA calibration.xlsx", sheet=1) #update path to data file location
# library(readxl)
# dat = read_xlsx("C:/Users/vtobias/Documents/Loach_Study_2023/Data_Raw/USFWS_Loach removal_eDNA calibration.xlsx")
dat = dryad_files_download('xxx')
# CalibrationExperiment_qPCR_RAW.csv

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
saveRDS(m_sc, "./Model_output/small_channel.rds") # export model results
# small_channel.rds provided in the Model_output folder on github is from the
#   run used in the manuscript. Each run will produce slightly different results.
m_sc = readRDS("./Model_output/small_channel.rds") # import model results (if already exported)

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
ggsave("./Figures/sc_model_plot2.jpeg",
       plot=sc_model_plot, 
       width=7, 
       height= 5.25) #1.5*aspect_ratio) # update directory as needed


# p-detect in small channel
# create a dataset with labels and predictor variables:
small_canal_pod <- data.frame(Location = rep("Small Channel", 10),
                              Target = c(rep("Ballyhoo", 5), rep("Anchovy", 5)),
                              Distance = c(500, 250, 100, 50, 10, 10, 500, 250, 100, 50),
                              dist_lc_m50  = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
                              dist_lc_m100 = c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0),
                              dist_lc_m250 = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
                              dist_lc_m500 = c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0),
                              TargetBHOO   = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                              n_techrep = rep(3, 10),
                              PoD_mean = NA, # empty column for storing predictions
                              PoD_low = NA,  # empty column for storing predictions
                              PoD_high = NA  # empty column for storing predictions
)

# Use the dataset above to predict the probability of detection
#  Cycle through the rows to use the data, then write the results to the
#  appropriate columns.
for(i in 1:nrow(small_canal_pod)){
  tmp <- est_p_detect(variable_levels = c(dist_lc_m50 = small_canal_pod$dist_lc_m50[i],
                                 dist_lc_m100= small_canal_pod$dist_lc_m100[i],
                                 dist_lc_m250= small_canal_pod$dist_lc_m250[i],
                                 dist_lc_m500= small_canal_pod$dist_lc_m500[i],
                                 TargetBHOO = small_canal_pod$TargetBHOO[i]),
             std_curve_alpha = 18.475, # for loach, from PL_MM assay standard curve
             std_curve_beta = -1.545,  # for loach, from PL_MM assay standard curve
             model_fit = m_sc, 
             n_rep = 3) # 3 reps = 3 filters
  small_canal_pod$PoD_mean[i] <- mean(tmp)
  small_canal_pod[i, c("PoD_low", "PoD_high")] <- quantile(tmp, probs = c(0.025, 0.975))
}


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
ggsave("./Figures/sc_pod_plot.jpeg",
       plot=sc_pod_plot, 
       width=5, 
       height= 7) #1.25*aspect_ratio)  # update directory

# PoD plot with error bars
facet_names = as_labeller(
  c(`Anchovy` = "Appx. 16g", 
    `Ballyhoo` = "Appx. 100g"))

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
        plot.subtitle = element_text(hjust=0.5), 
        legend.position = "bottom") +
  facet_wrap(~Target, 
             ncol= 2, 
             labeller=facet_names)


ggsave("./Figures/sc_pod_plot_errorbars.jpeg",
       plot=sc_pod_plot_eb, 
       width=10, 
       height= 6) #2*aspect_ratio) # update directory

### END ###
