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
source("./R_code/Field_Data_prep.R")

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
abline(v = 1.5, col = "grey")
text(2, -0.05, "Main Study", col = "grey", cex = 0.6)
text(0.8, -0.05, "Pilot Study", col = "grey", cex = 0.6)
axis(side = 1, 
     at = 1:6, 
     labels = c("Pilot", "Week 2", "Week 4", "Week 6", "Week 7", "Week 8"),
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
axis(side = 1, at = 1:6, 
     labels = c("Pilot", "Week 2", "Week 4", "Week 6", "Week 7", "Week 8"), 
     cex.axis = 0.75)
dev.off()

# extract values from the boxplot:
# rows are: lower whisker, 
#           the lower hinge, 
#           the median, 
#           the upper hinge and 
#           the extreme of the upper whisker
boxplot(total_length ~ round.f, 
        data = fish)$stats[3,]

# Linear model for differences in length among rounds
summary(lm(total_length ~ round.f,
           data = fish))

library(lubridate)
fish$date_set <- mdy(fish$date_set)
fish$doy <- yday(fish$date_set)
# linear model - length as a function of day - pooled both studies
summary(lm(total_length ~ doy,
           data = fish))
plot(fish$doy[-which(is.na(fish$total_length))],
     resid(lm(total_length ~ doy,
           data = fish)))
abline(h = 0, lty = 2, col = "grey")


plot(fish$doy,
       fish$total_length)
lines(min(fish$doy):max(fish$doy),
      predict(lm(total_length ~ doy,
                data = fish),
             newdata = data.frame(doy = min(fish$doy):max(fish$doy))))


# 5. Pilot Study Confusion Matrix ####
confusion_matrix(confDat$detectTrap,  # prediction vector
                 confDat$detectEDNA,  # "truth" vector
                 positive = "Yes",
                 # alpha = 0.95,
                 boot = FALSE)

# 6. Summarize zero catch ####

# probability of catching zero loach

# Pilot
predict(
  cpueP_tp,
  newdata = newdataPilot,
  type = "zprob"
)
# [1] 0.8615384

statPilotZprob <- function(df, inds) {
  predict(
    glmmTMB(loach_captured ~ 1,
            offset = log(total_time_fished),
            data = df[inds,],
            ziformula = ~1,
            family = truncated_poisson),
    newdata = newdataPilot,
    type = "zprob"
  )
}
set.seed(2018)
resPilotZprob <- boot(datPilot, 
                 statPilotZprob, 
                 R = 100)
CIPilotZprob <- setNames(as.data.frame(t(sapply(1:nrow(newdataPilot), 
                                           function(row) boot.ci(resPilotZprob, 
                                                                 conf = 0.95, 
                                                                 type = "basic", 
                                                                 index = row)$basic[, 4:5]))),
                    c("lower", "upper"))
predPilotZprob <- cbind.data.frame(newdataPilot, 
                              response = predict(
                                cpueP_tp,
                                newdata = newdataPilot,
                                type = "zprob"
                              ), 
                              CIPilotZprob)

# Main
predict(
  cpueM_tnb,
  newdata = newdataMain,
  type = "zprob"
)
# [1] 0.1168825 0.1168825 0.1168825 0.1168825 0.1168825

# Main study model with just one trap
# zprob <- data.frame(trapcount = 1:10,
#                     zprob = NA)
# for(i in 1:10){
#   cpueM_tnb_one <- glmmTMB(loach_captured ~ round.f, #trap round as a factor
#                      offset = log(total_time_fished),
#                      data = traps[which(traps$study == "Main" &
#                                           traps$trap_num %in% 1:i),],
#                      ziformula = ~ 1,
#                      family = truncated_nbinom2)
# 
# zprob[i, 2] <- predict(cpueM_tnb_one,
#         newdata = data.frame(round.f = newdataMain$round.f,
#                              total_time_fished = newdataMain$timeFishedMinTotal),
#         type = "zprob")[1]
# }
# 
# plot(zprob,
#      pch = 16, type = "b",
#      xlab = "Traps Included")

# bootstrap estimates of zprob using total counts of specified numbers of traps
# create dataset to contain results
zprob <- data.frame(trapcount = 1:10,
                    zprob = NA,
                    lcl = NA,
                    ucl = NA)
# specify statistics for the boot function to use
statZprob <- function(df, inds) {
  predict(
    glmmTMB(loach_captured ~ round.f, #trap round as a factor
                    offset = log(total_time_fished),
                    data = df[inds,],
                    ziformula = ~ 1,
                    family = truncated_nbinom2),
    newdata = data.frame(round.f = newdataMain$round.f,
                         total_time_fished = newdataMain$timeFishedMinTotal),
    type = "zprob"
  )
}

for(i in 1:10){
  print(i) # helps keep track of progress because this take a little while to run
  set.seed(2018)
  # create a new dataset like trapsTotals, but only including the specified traps
  df <- traps[which(traps$study == "Main" &
                      traps$trap_num %in% c(1:i)),] %>% 
    group_by(site_name) %>% 
    summarise(loach_captured = sum(loach_captured, na.rm = TRUE),
              trapCount = n(),
              total_time_fished = sum(total_time_fished, na.rm = TRUE),
              round.f = unique(round))
  # do the bootstrapping
  resZprob <- boot(df, 
                   statZprob, 
                   R = 100) 
  # calculate confidence intervals
  CIZprob <- setNames(as.data.frame(t(sapply(1:nrow(data.frame(round.f = newdataMain$round.f,
                                                               total_time_fished = newdataMain$timeFishedMinTotal)), 
                                             function(row) boot.ci(resZprob, 
                                                                   conf = 0.95, 
                                                                   type = "basic", 
                                                                   index = row)$basic[, 4:5]))),
                      c("lower", "upper"))
  # organize the results
  predZprob <- cbind.data.frame(trapCount = i, 
                                response = predict(
                                  glmmTMB(loach_captured ~ round.f, #trap round as a factor
                                          offset = log(total_time_fished),
                                          data = df,
                                          ziformula = ~ 1,
                                          family = truncated_nbinom2),
                                  newdata = data.frame(round.f = newdataMain$round.f,
                                                       total_time_fished = newdataMain$timeFishedMinTotal),
                                  type = "zprob"
                                ), 
                                CIZprob)
  # put the results into the container dataset
  zprob[i,] <- predZprob[1,] # could just do one round in the newData because they're all the same in this model
}

plotCI(x = zprob$trapcount,
       y = zprob$zprob,
       ui = zprob$ucl,
       li = zprob$lcl,
     pch = 16, 
     ylim = c(0, 1),
     xlab = "Traps included \n (Trap 1 to trap number indicated)",
     ylab = "Probability of catching zero loaches at a site")
plotCI(x = 1,
       y = predPilotZprob$response,
       ui = predPilotZprob$upper,
       li = predPilotZprob$lower,
       pch = 16,
       col = "red",
       ylim = c(0, 1),
       add = TRUE)
legend("topright",
       pch = 16,
       col = c("red", "black"),
       legend = c("Pilot", "Main"))



# checking whether we need separate zero models for the main study
# Answer is no.
summary(glmmTMB(loachTotal ~ round.f, #trap round as a factor
        offset = log(timeFishedMinTotal),
        data = datMain,
        ziformula = ~  round.f,
        family = truncated_nbinom2))
