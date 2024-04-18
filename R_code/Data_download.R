
#### 1. Information ####
# Project: Invasive Loach eDNA

# Purpose: 
#    This script downloads data from the Dryad data repository.

# Author:
#    Vanessa Tobias <vanessa_tobias@fws.gov>

# Date updated:
#    2024-04-18


#### 2. Setup ####
library(rdryad)
# Note that the developer for {rdryad} plans for it to be superseded by 
#     {deposits} in the future, but at the time this script was developed
#     the transition had not happened yet.


#### 3. Data ####

# 3.1 Pilot Study Data
datPilot <- dryad_files_download('')

# 3.2 eDNA Follow-Up Study Data
