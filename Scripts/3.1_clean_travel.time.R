## -------------------------------------------------------------------------
##
## Script name: 3.1_clean_travel.time.R
##
## Purpose of script: Clean the output from the travel time script that was run on a
##                    HPC computer.
##
## Author: Masumi Stadler
##
## Date Finalized: 2024-02-09
##
## Copyright (c) Masumi Stadler, 2026
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes:
##
##
## -------------------------------------------------------------------------

## Use R project with regular scripts, all paths are relative 

# Server set-up -----------------------------------------------------------
## Working directory is set from where the job is submitted
## Load library path, if on a server
# .libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )

# R-setup -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
pckgs <- list('plyr', 'tidyverse','data.table' # wrangling & programming
             ) # change as needed

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

## Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

## Load custom functions --------------------------------------------------
funs <- list.files("./Functions", full.names = T)
invisible(lapply(funs, source))

## Other set-up -----------------------------------------------------------
options(scipen = 6, digits = 4) # view outputs in non-scientific notation
 
## Parallel environment ---------------------------------------------------
## Server version
# cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
# registerDoMC(cores = cores)

## Personal version
# detectCores(); registerDoMC(); getDoParWorkers()
# numCores <- detectCores()
# cl <- makeCluster(numCores, type = "FORK")
 
# Read in data -----------------------------------------------------------

# Read in calculated travel time -----------------------------------------
y <- c(2015:2016)
mon <- c("6","8")

# create empty data frames to fill
out <- data.frame()

# Read in files and merge them together
for(i in 1:length(y)){
  for(j in 1:length(mon)){
  # Read in outputs
  # flow-weighted travel time aka water age
  cumtime <- readRDS(paste0("./Data/Traveltime/cumtime_",y[i], "-", mon[j],"_2024-02-09.rds")) %>% setDT()
  # stream network
  streamnet <-
    readRDS(paste0("./Objects/flacc3000_streamnet_wWRT_", y[i], "-", mon[j], "_final.rds"))
  
  m <- streamnet
  m <- m[!duplicated(m$pointid), ]
  m <- m[water.body == "Fluvial", flag_main.chan := 1]
  
  temp <- merge(m, cumtime, by = "pointid")
  
  temp[, water.body := factor(water.body, levels = c("Fluvial", "Reservoir", "Lake"))]
  temp[, travel.time_d := travel.time * 0.00069444] # convert from mins to day
  
  temp[, coord_x := as.numeric(coord_x)]
  temp[, coord_y := as.numeric(coord_y)]
  
  out <- rbind(out, temp)
  
}}


# # 2017-2018 --------------------------------------------------------------------------------------------------------
# y <- c(2017:2018)
# mon <- c("6","8","10")
# 
# # Read in files and merge them together
# for(i in 1:length(y)){
#   for(j in 1:length(mon)){
#     # Read in outputs
#     # flow-weighted travel time aka water age
#     cumtime <- readRDS(paste0("./Data/Traveltime/cumtime_",y[i], "-", mon[j],"_2024-02-09.rds")) %>% setDT()
#     # stream network
#     streamnet <-
#       readRDS(paste0("./Objects/flacc3000_streamnet_wWRT_", y[i], "-", mon[j], "_final.rds"))
#     
#     m <- streamnet
#     m <- m[!duplicated(m$pointid), ]
#     m <- m[water.body == "Fluvial", flag_main.chan := 1]
#     
#     temp <- merge(m, cumtime, by = "pointid")
#     
#     temp[, water.body := factor(water.body, levels = c("Fluvial", "Reservoir", "Lake"))]
#     temp[, travel.time_d := travel.time * 0.00069444] # convert from mins to day
#     
#     temp[, coord_x := as.numeric(coord_x)]
#     temp[, coord_y := as.numeric(coord_y)]
#     
#     out <- rbind(out, temp)
# }}
# 
# max.trav <- out[flag_mouth == 1,]
# 
# max.trav[, mouth.trav.time_d := travel.time * 0.00069444]
# write.table(max.trav, "./Data/Traveltime/travel.time_byyear_atmouth.csv", sep = ",", dec = ".", row.names = F)

# Assign values to samples ---------------------------------------------------------------------
# read in samples and their pointids
rest <- read.csv("./Data/Traveltime/restsamples_snapped_mainchan_2018_pointid.txt", sep = ";", dec = ",") %>%
  dplyr::select(sample.name = sample_nam, pointid) %>% setDT()

main <- read.csv("./Data/Traveltime/mainsamples_snapped_mainchan_2018_pointid.txt", sep = ";", dec = ",") %>%
  dplyr::select(sample.name = sample_nam, pointid) %>% setDT()

# # Read in samples from 2017, 2018 and their snapped pointids
# s1718 <- read.csv("./Data/Traveltime/20172018_samples_DOMMic_near3000_pointid.txt", sep = ";", dec = ",") %>%
#   dplyr::select(sample.name = sample_nam, pointid = flacc3000_) %>% setDT()

# gather the rest of meta data
mother.main <- read.csv("./Data/Summary/mother_mainsamples.csv", sep = ",", dec = ".") %>% setDT()
mother.rest <- read.csv("./Data/Summary/mother_restsamples.csv", sep = ",", dec = ".") %>% setDT()
mother.meta <- rbind(mother.main, mother.rest)

# combine the samples on main channel and the rest snapped to tributaries
meta <- rbind(main, rest)

# # add 2017/18 data into meta
# meta <- meta[s1718, pointid := i.pointid, on = .(sample.name)]

# correct a few manually
meta[pointid == 525799, pointid := 526769]
meta[pointid == 457095, pointid := 452871]

# merge the two
meta <- mother.meta[meta, , on = .(sample.name)]
# rm(main, rest, mother.main, mother.rest, mother.meta)
meta <- meta[!is.na(year),]
meta[, month := as.numeric(as.character(factor(campaign, levels = c(1,2,3), labels = c(6,8,10))))]

meta[sample.name %in% meta[sample.type.year == "Marine",]$sample.name, 
     pointid := 547615] # overwrite marine samples with mouth

nrow(meta) #1761

meta[, Year := as.factor(year)]
meta[, Month := as.factor(month)]
out[, Year := as.factor(Year)]
out[, Month := as.factor(Month)]


meta[out, c("flacc_km2", "strahler.order", "dis_m3s", "vel_ms",
                                    "wrt_min", "travel.time_min", "travel.time_days") :=
                         list(i.flacc_km2, i.strah_ord, i.dis_m3s, i.vel_ms.wb, 
                              i.wrt_min.wb, i.travel.time, i.travel.time_d), on = c("Year", "Month", "pointid")]


meta$sample.type.year <- factor(meta$sample.type.year, levels = c("Soil","Sediment",
                                                                  "Soilwater","Hyporheic",
                                                                  "Wellwater","Stream", "Tributary",
                                                                  "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                  "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                  "Marine"),
                                labels = c("Soil","Sediment",
                                           "Soilwater","Soilwater",
                                           "Groundwater","Stream", "Tributary",
                                           "Riverine Lakes", "Headwater Ponds", "Lake", "Lake",
                                           "Upriver",# "RO3", "RO2", "RO1","Deep",
                                           "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
                                           "Downriver",
                                           "Estuary"))

# save meta file
write.table(meta, paste0("./Data/Summary/meta_file_traveltime_", Sys.Date(), ".csv"), sep = ",", dec = ".", row.names = F)





sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: x86_64-redhat-linux-gnu
# Running under: Fedora Linux 43 (Workstation Edition)
# 
# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.12.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: America/New_York
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] doParallel_1.0.17   iterators_1.0.14    foreach_1.5.2       ggpubr_0.6.3       
# [5] plyr_1.8.9          lubridate_1.9.5     forcats_1.0.1       stringr_1.6.0      
# [9] dplyr_1.2.0         purrr_1.2.1         readr_2.2.0         tidyr_1.3.2        
# [13] tibble_3.3.1        ggplot2_4.0.2       tidyverse_2.0.0     data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] generics_0.1.4     rstatix_0.7.3      stringi_1.8.7      hms_1.1.4         
# [5] magrittr_2.0.4     grid_4.5.2         timechange_0.4.0   RColorBrewer_1.1-3
# [9] backports_1.5.0    Formula_1.2-5      scales_1.4.0       CoprManager_0.5.8 
# [13] codetools_0.2-20   abind_1.4-8        cli_3.6.5          rlang_1.1.7       
# [17] withr_3.0.2        tools_4.5.2        tzdb_0.5.0         ggsignif_0.6.4    
# [21] broom_1.0.12       vctrs_0.7.1        R6_2.6.1           lifecycle_1.0.5   
# [25] car_3.1-5          pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6      
# [29] glue_1.8.0         Rcpp_1.1.1         tidyselect_1.2.1   rstudioapi_0.18.0 
# [33] farver_2.1.2       carData_3.0-6      labeling_0.4.3     compiler_4.5.2    
# [37] S7_0.2.1     