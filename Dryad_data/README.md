# Data from: Leveraging Environmental DNA (eDNA) to Optimize Invasive Fish Eradication Efforts

\[[https://doi.org/10.5061/dryad.n8pk0p342](https://doi.org/10.5061/dryad.n8pk0p342)]

## General information

### Author Information

Corresponding author:
Jennie J. Wiggins
Fish Biologist, Lodi Fish and Wildlife Office
US Fish and Wildlife Service
[jennie_wiggins@fws.gov](mailto:jennie_wiggins@fws.gov)

Co-Authors:
Vanessa Tobias
Mathematical Statistician, Lodi Fish and Wildlife Office
US Fish and Wildlife Service
[vanessa_tobias@fws.gov](mailto:vanessa_tobias@fws.gov)

Erika F. Holcombe
Supervisory Fish Biologist, Lodi Fish and Wildlife Office
US Fish and Wildlife Service
[erika_holcombe@fws.gov](mailto:erika_holcombe@fws.gov)

Katie Karpenko
Research Associate, Genidaqs Laboratory
Cramer Fish Sciences
[katie.karpenko@fishsciences.net](mailto:katie.karpenko@fishsciences.net)

Eric R. Huber
Supervisory Fish Biologist, Lodi Fish and Wildlife Office
US Fish and Wildlife Service
[eric_huber@fws.gov](mailto:eric_huber@fws.gov)

Andrew C. Goodman
Fish Biologist, Lodi Fish and Wildlife Office
US Fish and Wildlife Service
[andrew_goodman@fws.gov](mailto:andrew_goodman@fws.gov)

Date of data collection: Pilot Study- March 2022; Main Study (including calibration experiment)- February 2023 – June 2023

Geographic location of data collection: All data were collected at San Luis National Wildlife Refuge (Merced County, California, USA)

## Data and file overview

This dataset contains 8 CSV files:
•	All_Variable_descriptions.csv
•	catch_loach_data.csv
•	loach_length_data.csv
•	eDNA_detections_site_data.csv
•	Resampled_eDNA_detections_site_data.csv
•	PilotStudy_qPCR_RAW.csv
•	CalibrationExperiment_qPCR_RAW.csv
•	MainStudy_qPCR_RAW.csv

### File descriptions:

All_Variable_descriptions.csv contains 4 columns and 139 rows of descriptions for all the variables used within the 7 other CSV files associated with this dataset.

catch_loach_data.csv contains 21 variables and 897 (898 with column headers) rows of data related to loach capture per minnow trap collected during trapping efforts for the Pilot Study and the Main Study. This file only contains data related to trapped sites.Columns are study, round, site_name, trap_name, habitat_type, trap_num, latitude, longitude, date_set, time_set, date_pulled, time_pulled, total_time_fished, loach_captured, dissolved_oxygen, conductivity, temperature, turbidity, depth, velocity, notes. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file. Environmental data has not been scrubbed for outliers.

loach_length_data.csv contains 16 variables and 676 (677 with column headers) rows of length data for each loach collected during trapping efforts for the Pilot Study and the Main Study. Columns are study, round, site_name, trap_name, habitat_type, trap_num, latitude, longitude, date_set, time_set, date_pulled, time_pulled, species_name, species_code, total_length, notes. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file.

eDNA_detections_site_data.csv contains 16 variables and 377 (378 with column headers) rows of categorical and environmental data related to all sites sampled for eDNA and/or trapped during the Pilot Study or Main Study. Columns are study, round, site_name, habitat_type, latitude, longitude, date, sampled_eDNA, detection, trapped, dissolved_oxygen, conductivity, temperature, turbidity, depth, velocity. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file. Environmental data has not been scrubbed for outliers.

Resampled_eDNA_detections_site_data.csv contains 15 variables and 38 (39 with column headers) rows of categorical and environmental data related to sites that were re-sampled for eDNA if loach were captured there during trapping efforts for the Main Study. Columns are study, round, site_name, habitat_type, latitude, longitude, date, detection, dissolved_oxygen, conductivity, temperature, turbidity, depth, velocity, notes. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file. Environmental data has not been scrubbed for outliers.

PilotStudy_qPCR_RAW.csv contains 20 variables and 540 (541 with column headers) rows of qPCR results for all filters collected during the Pilot Study. Columns are projectID, project_type, extraction_order, date, site_name, filter_replicate, FilterID, vol_filtered_ml, location, location_broad, extraction_date, control, month, year, field_notes, technical_replicate, Cq, target, relative_conc_ng/uL, relative_conc_weighted_by_vol. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file.

CalibrationExperiment_qPCR_RAW.csv contains 33 variables and 630 (631 with column headers) rows of qPCR results for all filters collected during the eDNA calibration experiment of the Main Study. Columns are projectID, project_type, extraction_order, date, name, filter_replicate, FilterID, vol_filtered_ml, location, location_broad, extraction_date, lat, lon, control, month, year, field_notes,dist_lc_m, lc_deploy_time, lc_remove_time, lc_spp_1, lc_spp_2, spp_1_wt_g, spp_2_wt_g , vel_m_s, vel_ft_sec, turb_ntu, depth_m,  technical_replicate, Cq, relative_conc_ng/uL, target. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file.

MainStudy_qPCR_RAW.csv contains 22 variables and 2674 (2675 with column headers) rows of qPCR results for all filters collected during the Main Study. Columns are projectID, project_type, extraction_order, date, site_name, filter_replicate, FilterID, vol_filtered_ml, location, location_broad, extraction_date, lat, lon, control, month, year, field_notes, technical_replicate, Cq, target, relative_conc_ng/uL, relative_conc_weighted_by_vol. A complete description of all the variables used in this file can be found in the All_Variable_descriptions.csv file.

### Methodological information

Please refer to the associated paper for a more detailed overview of data collection, processing, and analyses for this dataset.

eDNA sample collection: Environmental DNA samples were collected by simultaneously filtering water through three Millipore® Sterivex™-HV 0.45 μm sterile PVDF membrane filter units (Millipore Sigma, Burlington, MA, USA) using three Masterflex® L/S Easy-Load II peristaltic pumps. Environmental DNA field controls were collected at the start of each sampling day to confirm that field equipment was free of contamination.

Pilot Study: We selected 127 random sampling sites made up of an approximately equal number of habitat types termed ponds (n=63) and canals (n=64). One baited minnow trap was set at each of the selected locations for approximately 24 h. A random subset of sites were also sampled for eDNA, including 19 canal sites and nine pond sites. We measured the total length (TL) of all trapped loaches. Species identification was genetically confirmed for each specimen captured resulting in 17 *Paramisgurnus dabryanus* (Large-Scale Loach) and one *Misgurnus mizolepis* (Fine-Scale Loach).

Calibration experiment: To better understand eDNA transport and persistence within the sampling area, we conducted a fixed-point DNA experiment in a small shallow canal and a large deep canal. We created a fixed source of eDNA using known quantities of commercially available frozen *Engraulis mordax* (Northern Anchovy; 16 g) and *Hemiramphus brasiliensis* (Ballyhoo; 100 g). Velocity and turbidity were measured prior to sampling using a Global Water FP111 Flow Probe (YSI Inc., Yellow Springs, OH, USA) and a Hach 2100Q Turbidity Meter (Hach Co., Loveland, CO, USA), respectively. After the calculated dispersal time had passed, we collected eDNA samples simultaneously in triplicate at six set distances from the tethered fish. Once complete, eDNA sources were removed from each canal. The eDNA collection process was repeated at the same locations the following day to determine if eDNA from the proxy species remained in the system after 24 h.

Main Study: We selected 300 sampling sites comprised of only canal habitat. Environmental DNA sampling and loach trapping alternated weekly to permit strategic deployment of traps based on the previous week’s eDNA results; areas containing the highest concentrations of loach eDNA received priority for trapping. Two additional weeks of trapping only were conducted after the eDNA filter budget was exhausted. At sites where loach eDNA was detected, we deployed 10 minnow traps over a total of 260 m moving upstream. Traps soaked for approximately 24 h before they were removed. Sites where loach were captured, were resampled for eDNA the following week (designated by an ‘RX’ in the site name). We measured the TL of all trapped loaches. Species identification for the loaches collected during this study have not been genetically confirmed and thus have been assigned a generic species name and code at the time of this publication. Environmental data were collected at each site one time before eDNA sampling and one time before trapping.

Code/Software
The R script files for reproducing the analyses in the associated manuscript can be found at [https://github.com/USFWS/LFWO-Loach-eDNA.git](https://github.com/USFWS/LFWO-Loach-eDNA.git). The scripts were created using R version 4.3.1 "Beagle Scouts".
