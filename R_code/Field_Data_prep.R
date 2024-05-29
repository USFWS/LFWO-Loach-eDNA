
#### 1. Information ####
# Project: Invasive Loach eDNA

# Purpose: 
#   This script downloads data from the Dryad data repository associated with 
#   Wiggins et al. Leveraging Environmental DNA (eDNA) to Optimize Invasive Fish 
#   Eradication Efforts. We provide an additional scripts that 
#    1. analyze the data downloaded by this script to produce figures and tables
#       related to loach catch and size in the manuscript and
#    2. calibrate the eDNA detection data used in this script and that generates 
#       the graphs associated with eDNA calibration in the manuscript.

# Author:
#    Vanessa Tobias <vanessa_tobias@fws.gov>

# Date updated:
#    2024-05-17


#### 2. Setup ####
library(rdryad)
# Note that the author of this script is aware that the developer for {rdryad} 
#     plans for it to be superseded by {deposits} in the future, but at the 
#     time this script was developed the transition had not happened yet.
#     Potential future users of this script may need to update the code once
#     this transition has happened.


#### 3. Download Data ####

# Datasets are available on Dryad:
# https://doi.org/10.5061/dryad.n8pk0p342
# use this to get Dryad's ID numbers:
dryad_dataset("10.5061/dryad.n8pk0p342")
# use this to get version numbers:
dryad_dataset_versions("10.5061/dryad.n8pk0p342")

##### 3.1 eDNA ####

# The code below will download data directly from Dryad once the data is made
#      public. Peer reviewers can use the link provided by Dryad to download 
#      the data. Filenames are provided below.

# Site-level eDNA information 
edna <- dryad_files_download('xxx')
# eDNA_detections_site_data.csv

## 3.2 Traps ####
traps <- dryad_files_download('xxx')
# catch_loach_data.csv

## 3.3 Fish ####
fish <- dryad_files_download('xxx')
# loach_length_data.csv


#### 3. Prepare Data ####

##### 3.1 eDNA ####



## 3.2 Traps ####
# calculate CPUE as loaches/minute
traps$cpue_Min <- traps$loach_captured/traps$total_time_fished

# change round labels so they're consistent
traps$round[which(is.na(traps$round))] <- "Pilot"
traps$round[which(traps$round == "1")] <- "Round1"
traps$round[which(traps$round == "2")] <- "Round2"
traps$round[which(traps$round == "3")] <- "Round3"
traps$round[which(traps$round == "Extra_1")] <- "Extra1"
traps$round[which(traps$round == "Extra_2")] <- "Extra2"

# make a factor variable for the trap rounds
traps$round.f <- factor(traps$round,
                        levels = c("Pilot",
                                   "Round1", 
                                   "Round2", 
                                   "Round3", 
                                   "Extra1",
                                   "Extra2"))

# Create a summary dataset for all traps at each site
trapsTotals <- traps %>% 
  group_by(site_name) %>% 
  summarise(loachTotal = sum(loach_captured, na.rm = TRUE),
            trapCount = n(),
            timeFishedMinTotal = sum(total_time_fished, na.rm = TRUE),
            round = unique(round))

# make a factor variable for the trap rounds
trapsTotals$round.f <- factor(trapsTotals$round,
                              levels = c(
                                "Pilot",
                                "Round1", 
                                "Round2", 
                                "Round3", 
                                "Extra1",
                                "Extra2"))

confDat <- merge(edna,
                 trapsTotals[which(trapsTotals$round == "Pilot"),
                             c("site_name", "loachTotal")],
                 by = "site_name",
                 all.y = TRUE)
# remove rows where no eDNA samples were taken
confDat <- confDat[which(confDat$sampled_eDNA == "Yes"),]

confDat$detectTrap <- "Yes"
confDat$detectTrap[which(confDat$loachTotal == 0)] <- "No"

confDat <- data.frame(detectTrap = as.factor(confDat$detectTrap),  # prediction vector
                      detectEDNA = as.factor(confDat$detection)) #eDNA detection



## 3.3 Fish ####

# change round labels so they work as labels on graphs
fish$round[which(is.na(fish$round))] <- "Pilot"
fish$round[which(fish$round == "1")] <- "Round1"
fish$round[which(fish$round == "2")] <- "Round2"
fish$round[which(fish$round == "3")] <- "Round3"
fish$round[which(fish$round == "Extra_1")] <- "Extra1"
fish$round[which(fish$round == "Extra_2")] <- "Extra2"

fish$round.f <- factor(fish$round,
                       levels = c("Pilot",
                                  "Round1", 
                                  "Round2", 
                                  "Round3", 
                                  "Extra1",
                                  "Extra2"))


