# Purpose:
#   This script reproduces the analysis and graphs in Wiggins et al. 
#   Leveraging Environmental DNA (eDNA) to Optimize Invasive Fish 
#   Eradication Efforts. We provide an additional script that calibrates the 
#   eDNA detection data used in this script and that generates the graphs 
#   associated with eDNA calibration in the manuscript.

# Code author: 
#   Vanessa Tobias (USFWS; vanessa_tobias@fws.gov)

# Last updated: 
#   2024-05-17

# 1. Setup ####
# Load Packages
library(tidyverse) # data summary functions
library(glmmTMB)   # zero-inflated models for CPUE
library(DHARMa)    # model checking on glmmTMB models
library(boot)      # confidence intervals on CPUE models
library(plotrix)   # graphing
library(qwraps2)   # confusion matrix

# 2. Data ####

# Get data
source("./R_code/Data_download.R")

# 3. CPUE Models ####

## 3.1 Pilot Study ####
datPilot <- traps[which(traps$study == "Pilot" &
                          traps$habitat_type == "canal"),]

# Select best distribution
cpueP_p <- glmmTMB(loach_captured ~ 1,
                   offset = log(total_time_fished),
                   data = datPilot,
                   family = poisson)

cpueP_nb <- glmmTMB(loach_captured ~ 1,
                    offset = log(total_time_fished),
                    data = datPilot,
                    family = nbinom2())

cpueP_tp <- glmmTMB(loach_captured ~ 1,
                    offset = log(total_time_fished),
                    data = datPilot,
                    ziformula = ~1,
                    family = truncated_poisson)

cpueP_tnb <- glmmTMB(loach_captured ~ 1,
                     data = datPilot,
                     ziformula = ~1,
                     offset = log(total_time_fished),
                     family = truncated_nbinom2)

AIC(cpueP_p, cpueP_nb, cpueP_tp, cpueP_tnb)
# best = tp

#model checking
plot(DHARMa::simulateResiduals(cpueP_tp))

# Bootstrap CIs for mean CPUE
newdataPilot = data.frame(total_time_fished = c(1440))

statPilot <- function(df, inds) {
  predict(
    glmmTMB(loach_captured ~ 1,
            offset = log(total_time_fished),
            data = df[inds,],
            ziformula = ~1,
            family = truncated_poisson),
    newdata = newdataPilot,
    type = "response"
  )
}
set.seed(2018)
resPilot <- boot(datPilot, 
                 statPilot, 
                 R = 100)
CIPilot <- setNames(as.data.frame(t(sapply(1:nrow(newdataPilot), 
                                           function(row) boot.ci(resPilot, 
                                                                 conf = 0.95, 
                                                                 type = "basic", 
                                                                 index = row)$basic[, 4:5]))),
                    c("lower", "upper"))
predPilot <- cbind.data.frame(newdataPilot, 
                              response = predict(
                                cpueP_tp,
                                newdata = newdataPilot,
                                type = "response"
                              ), 
                              CIPilot)
names(predPilot)[which(names(predPilot) == "total_time_fished")] <- "timeFishedMinTotal"
predPilot$round <- "Pilot"

## 3.2 Main Study ####
# subset the site-level data on catch totals to the main study
datMain <- trapsTotals[which(trapsTotals$round != "Pilot"),]

cpueM_p <- glmmTMB(loachTotal ~ round.f, #trap round as a factor
                   offset = log(timeFishedMinTotal),
                   data = datMain,
                   ziformula = ~ 1,
                   family = poisson)

cpueM_nb <- glmmTMB(loachTotal ~ round.f, #trap round as a factor
                    offset = log(timeFishedMinTotal),
                    data = datMain,
                    ziformula = ~ 1,
                    family = nbinom2)

cpueM_tp <- glmmTMB(loachTotal ~ round.f, #trap round as a factor
                    offset = log(timeFishedMinTotal),
                    data = datMain,
                    ziformula = ~ 1,
                    family = truncated_poisson)


cpueM_tnb <- glmmTMB(loachTotal ~ round.f, #trap round as a factor
                     offset = log(timeFishedMinTotal),
                     data = datMain,
                     ziformula = ~ 1,
                     family = truncated_nbinom2)
# AIC(cpueE_p, cpueE_nb, cpueE_tp, cpueE_tnb, cpueE_p_cr, cpueE_tnb_temp, cpueE_tnb_tr)

AIC(cpueM_p, cpueM_nb, cpueM_tp, cpueM_tnb)
# cpueM_tnb

newdataMain = data.frame(round.f = factor(c("Round1", 
                                            "Round2", 
                                            "Round3", 
                                            "Extra1",
                                            "Extra2"),
                                          levels = c("Round1", #levels have to be in the same order as the original dataset
                                                     "Round2", 
                                                     "Round3", 
                                                     "Extra1",
                                                     "Extra2")),
                         timeFishedMinTotal = c(1440, 1440, 1440, 1440, 1440))

statMain <- function(df, inds) {
  #model <- formula(repro ~ fertilizer + level | fertilizer * level)
  predict(
    glmmTMB(loachTotal ~ round.f, #trap round as a factor
            offset = log(timeFishedMinTotal),
            data = datMain[inds,],
            ziformula = ~ 1,
            family = truncated_nbinom2),
    newdata = newdataMain,
    type = "response"
  )
}
set.seed(2018)
resMain <- boot(datMain, 
                statMain, 
                R = 100)
CIMain <- setNames(as.data.frame(t(sapply(1:nrow(newdataMain), 
                                          function(row) boot.ci(resMain, 
                                                                conf = 0.95, 
                                                                type = "basic", 
                                                                index = row)$basic[, 4:5]))),
                   c("lower", "upper"))
predMain <- cbind.data.frame(newdataMain, 
                             response = predict(
                               cpueM_tnb,
                               newdata = newdataMain,
                               type = "response"
                             ), 
                             CIMain)

predMain$round <- as.character(predMain$round.f)


predDF <- rbind(predPilot,
                predMain[, -1])
predDF$round.f <- factor(predDF$round,
                         levels = c("Pilot",
                                    "Round1", 
                                    "Round2", 
                                    "Round3", 
                                    "Extra1",
                                    "Extra2"))

png("./Figures/Fig_7.png",
    height = 6, width = 6, units = "in", res = 300)
par(cex = 1.25)
plotCI(x = 1:6,
       y = predDF$response,
       ui = predDF$upper,
       li = predDF$lower,
       pch = 16,
       xlim = c(0.5, 6.5),
       ylim = c(0, 4),
       xlab = "",
       ylab = "Mean catch per trap per day",
       xaxt = "n")
axis(side = 1, 
     at = 1:6, 
     labels = predDF$round,
     cex.axis = 0.7)
dev.off()

# 4. Total Length ####

png("./Figures/Fig_6.png",
    width = 6, height = 4, units = "in", res = 300)
par(cex = 1.25,
    mar = c(3.1, 4.1, 1.1, 2.1))
boxplot(total_length ~ round.f, 
        data = fish,
        xlab = "",
        ylab = "Total Length (mm)",
        ylim = c(15, 190), 
        xaxt = "n")
abline(v = 1.5, col = "grey")
text(2.3, 18, "Main Study", col = "grey", cex = 0.6)
text(0.8, 18, "Pilot Study", col = "grey", cex = 0.6)
axis(side = 1, at = 1:6, labels = levels(fish$round.f), cex.axis = 0.75)
dev.off()

# Linear model for differences in length among rounds
summary(lm(total_length ~ round.f,
           data = fish))

# 5. Pilot Study Confusion Matrix ####
confusion_matrix(confDat$detectTrap,  # prediction vector
                 confDat$detectEDNA,  # "truth" vector
                 positive = "Yes",
                 # alpha = 0.95,
                 boot = FALSE)