## -------------------------------------------------------------------------
##
## Script name: 5_MOclassification.R
##
## Purpose of script: This script is the main script for bacterial analyses.
##                    We identify reactive versus unreactive operational taxonomic units
##                    using a decision tree.
##
##
## Author: Masumi Stadler
##
## Date Finalized: 2024-02-12
##
## Copyright (c) Masumi Stadler, 2026
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes:  At sections in this script, one moves into scripts
##         5.1 and 5.2 and then, their outputs are continued to be
##         analyzed in this script.
##         The script is written in this way, because some sections
##         were run on a high-performance computer provided by
##         Compute Canada.
##         This script also produces Table S2.
##
## -------------------------------------------------------------------------

## Use R project with regular scripts, all paths are relative 

# Server set-up -----------------------------------------------------------
## Working directory is set from where the job is submitted
## Load library path, if on a server
# .libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )

# R-setup -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
pckgs <- list(
  "phyloseq", "plyr", "tidyverse", "data.table", # wrangling
  "ggpubr",  "cowplot", "plotly", "gridExtra", # plotting, "ggnewscale","ggpmisc", 
  "kableExtra",# "xlsx", # making tables for manuscript (LaTeX) and export as excel (for ISME)
  "doMC", # parallel computing
  "vegan", "ape", "ade4", "rstatix", "mgcv") # statistics

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

# for Bioconductor packges run something like:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

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

# Personal version
detectCores(); registerDoMC(); getDoParWorkers()
numCores <- detectCores()
cl <- makeCluster(numCores, type = "FORK")

## Reproducibility set-up ------------------------------------------------
# Set seed for session and reproducibility of permutations
# (NMDS were done, but not part of main workflow anymore)
set.seed(3)

# Read in data -----------------------------------------------------------
# Output of first paper. Workflow in https://github.com/CarBBAS/Paper_Stadler-delGiorgio_ISMEJ_2021
# Output of script: 2_analysis.R

# ASV CSS transformed table
otu.tab <- select_newest("./Data/Microbial", "201520162017_fin_css_otu99_phantomcor_no.dups_")
otu.tab <- as.matrix(read.csv(
  otu.tab,
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))

# row orders need to match between tax.tab and otu.tab
otu.tab <- otu.tab[order(row.names(otu.tab)),]

# Output of first paper. Workflow in https://github.com/CarBBAS/Paper_Stadler-delGiorgio_ISMEJ_2021
# Output of script: 1_correctByBin_css.R
# Taxonomy table
tax.tab <- select_newest("./Data/Microbial", "201520162017_tax_otu99_table_no.dups_")
tax.tab <-
  as.matrix(read.csv(
    tax.tab,
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ))
# orders need to match between tax.tab and otu.tab
tax.tab <- tax.tab[order(row.names(tax.tab)),]

# Output of first paper. Workflow in https://github.com/CarBBAS/Paper_Stadler-delGiorgio_ISMEJ_2021
# Output of script: 1_correctByBin_css.R
# Meta data
met.df <-
  select_newest(path = "./Data/Microbial", file.pattern = "201520162017_meta_otu99_no.dups_")
met.df <-
  read.csv(
    met.df,
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ) %>% setDT(keep.rownames = "row.names")


## Sample renaming ---------------------------------------------------------------------------------------
# add travel time
meta.all <- read.csv("./Data/Summary/meta_file_traveltime_2024-02-10.csv", sep = ",", stringsAsFactors = F) %>% setDT()

# omit sediment samples
met.df <- met.df[sample.type.year != "Sediment",]
met.df <- met.df[year <= 2016,]
met.df <- met.df[campaign <= 2,]

# omit estuary
met.df <- met.df[sample.type.year != "Marine",]

meta.all[, mic.name := sample.name]

# rename name to match taxonomy meta file to travel time meta
meta.all[str_detect(mic.name, "^RO"), mic.name := str_replace(mic.name, "[-]", "")]
meta.all[str_detect(mic.name, "^RO"), mic.name := str_replace(mic.name, "[_]", ".")]
meta.all[str_detect(mic.name, "hypo"), mic.name := str_replace(mic.name, "hypo", ".Hypo")]
meta.all[str_detect(mic.name, "^T49"), mic.name := str_replace(mic.name, "T49", "TR49")]
meta.all[str_detect(mic.name, "^T5"), mic.name := str_replace(mic.name, "T5", "TR5")]
meta.all[str_detect(sample.name, "^L32"), mic.name := str_replace(mic.name, "[_]", "")]
meta.all[str_detect(sample.name, "^L33"), mic.name := str_replace(mic.name, "[.]", "")]
meta.all[sample.name == "RO2-111hypo", mic.name := "RO2111.90m"]
meta.all[sample.name == "T49_C2", mic.name := "TR49"]
meta.all[sample.name == "L330_C1", mic.name := "L330.C1"]
meta.all[sample.name == "L330_C2", mic.name := "L330"]
meta.all[sample.name == "LR09.1", mic.name := "LR09"]
meta.all[sample.name == "SWPR02.1", mic.name := "SWPR02"]

# which don't match?
met.df$dr_match_name[!(met.df$dr_match_name %in% meta.all$mic.name)]

setDT(met.df)
# merge travel time to mic meta data
met.df[meta.all, c("Season","pointid", "strahler.order","flacc_km2",
                   "vel_ms","dis_m3s","wrt_min", "travel.time_d") := # "velocity_ms", 
         list(i.Season, i.pointid, i.strahler.order, i.flacc_km2,
              i.vel_ms, i.dis_m3s,  i.wrt_min, i.travel.time_days), on = c("dr_match_name==mic.name")] #i.velocity_ms,


setDF(met.df, rownames = met.df$row.names) # convert back to data frame

## Define habitats -----------------------------------------------------------------------------------------
# merge some sample types
met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheic",
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                      "Marine","Unknown"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Soilwater",
                                             "Groundwater","Stream", "Tributary",
                                             "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
                                             "Upriver",
                                             "Reservoirs","Reservoirs", "Reservoirs","Reservoirs", "Downriver",
                                             "Estuary","Unknown"))

met.df$Season <- factor(met.df$Season, levels = c("Spring","Summer","Autumn"),
                        labels = c("Spring","Summer","Autumn"))

met.df$dna_type <- factor(met.df$dna_type, levels = c("DNA","cDNA"),
                          labels = c("DNA","RNA"))

# save coordinates
out <- met.df %>% dplyr::select(seq_name, dr_match_name, dna_type, season, year, sample.type.year, lat, long)
write.table(out, "./Output/coordinates_MO.csv", sep = ",", row.names = F)

# Table S2 ----------------------------------------------------------------
# Number of samples used to create spatial models of microbial OTUs and DOM MF per campaign
# make sure we only count those samples we are using
temp <- met.df[met.df$row.names %in% row.names(otu.tab),] %>% setDT()

# get microbial N
temp <- temp[,.(n = .N), by = .(dna_type, year, season)] %>% 
  arrange(dna_type, year, season) %>%
  dplyr::rename(dataset = dna_type) %>%
  filter(season != "Autumn")

# load DOM N
temp.d <- read.csv("./Data/Summary/FT_samplesize.csv", sep = ",", stringsAsFactors = F) %>%
  mutate(dataset = "DOM") %>%
  dplyr::rename(season = Season) %>%
  filter(season != "Autumn")

# Merge
temp <- bind_rows(temp, temp.d) %>% setDT()
# format and re-order for table
temp <- dcast(temp, year + season ~ dataset, value.var = "n")
temp <- temp %>% dplyr::select(Year = year,
                       Season = season,
                       DOM, DNA, RNA)
temp[is.na(RNA),"RNA"] <- 0

# format output in LaTeX
knitr::kable(temp, "latex", booktabs = T) %>%
  kable_styling(position = "center", full_width = F)

rm(temp, temp.d)

temp <- met.df[met.df$row.names %in% row.names(otu.tab),] %>% 
  filter(season != "Autumn") %>%
  setDT()

# get microbial N
temp <- temp[,.(n = .N), by = .(dna_type, year, sample.type.year)] %>% 
  arrange(dna_type, year, sample.type.year) %>%
  dplyr::rename(dataset = dna_type)
  

## Converge all data and prepare for processing -------------------------------------------------------
# Construct phyloseq object
pb <- phyloseq(otu_table(otu.tab[row.names(otu.tab) %in% row.names(met.df),], taxa_are_rows = F),
               sample_data(met.df),
               tax_table(tax.tab))

pb <- prune_taxa(!taxa_sums(pb) == 0, pb)
pb <- prune_samples(!sample_sums(pb) == 0, pb)

# format into long
# melt community matrix for tidy data set
commat <- melt.data.table(
  setDT(as.data.frame(otu_mat(pb)), keep.rownames = "Sample"), #otu_mat(pb)
  #scale(otu_mat(pb), center = F, scale = T))
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

# Some stats
minhead(otu_mat(pb))
nrow(otu_mat(pb)) # 410 samples
ncol(otu_mat(pb)) # 10589 OTUs

# merge back with meta data
mic <- commat[setDT(met.df %>% dplyr::select(seq_name, dr_match_name, replicate, dna_type, seq_depth, year,
                                             Season, sample.type.year, sample.name, flacc_km2, dis_m3s,
                                             wrt_min,
                                             travel.time_d), keep.rownames = "Sample"),
              , on = .(Sample)]

## Remove outliers -------------------------------------------------------------------------------------------
# detect extreme group outliers, except zeros
mic[reads > 0, ext_outlier := is_extreme(reads), by = .(year, Season, dna_type, OTU)]

# sanity check
any(mic[is.na(ext_outlier),]$reads > 0)
# overwrite all zeros with FALSE
mic[reads == 0, ext_outlier := FALSE]

# check how many outliers
nrow(mic[ext_outlier == TRUE,]) * 100 / nrow(mic) # 0.53% outliers

# remove extreme outliers
clean.df <- mic[ext_outlier == FALSE,]

## Separate absent and present MF ---------------------------------------------------------------------------
# Define absence
clean.df[, ID := paste(year, Season, dna_type, OTU, sep = "_")] # year,
clean.df[, n := .N, by = .(ID)] # number of samples per category
clean.df[, n.obs := nrow(.SD[reads > 0,]), by = .(ID)] # number of samples with actual observation

#initiate absence column
clean.df[n.obs == 0, PA := "Absent", by = .(ID)]

#how many microbes that are absent?
length(levels(factor(clean.df[PA == "Absent",]$OTU))) # 94258
length(levels(factor(clean.df$OTU))) # 10589
# percentage lost
length(levels(factor(clean.df[PA == "Absent",]$OTU))) * 100 / length(levels(factor(clean.df$OTU)))
# 87.43% absent
# this doesn't mean that all OTUs never appear in this dataset. This means that only 11% appear in all years, campaigns
# and DNA vs RNA

# Overwrite the rest with present
clean.df[is.na(PA), PA := "Present"]
# extract only df with present OTUs
present <- clean.df[PA == "Present",]
levels(factor(present$sample.type.year))

## Define travel time for non estimated habitats ----------------------------------------------------------
present[sample.type.year == "Soil", travel.time_d := -50]
present[sample.type.year == "Soilwater", travel.time_d := -50]
present[sample.type.year == "Groundwater", travel.time_d := -50] #0.0000001

## Transform variables and remove zero observations --------------------------------------------------------
present[reads == 0, reads := NA]

# remove weird sample types
present <- present[sample.type.year != "Unknown", ]
present[, ID := paste(year, Season, dna_type, OTU, sep = ".")]

# arrange so that the data frame follows travel time
present <- present %>% group_by(ID) %>% dplyr::arrange(travel.time_d)

#remove samples with no travel time or season
setDT(present)
present <- present[!is.na(travel.time_d),]
present <- present[!is.na(Season),]

absent <- clean.df[PA == "Absent",]

# save
# saveRDS(clean.df, paste0("./Objects/allMO_long_",Sys.Date(),".rds"))
# saveRDS(present, paste0("./Objects/presentMO_long_", Sys.Date(), ".rds"))
# saveRDS(absent, paste0("./Objects/absentMO_long_", Sys.Date(),".rds"))

## Z-Scale -------------------------------------------------------------------------------------------------
temp <- dcast(present, Sample ~ OTU, value.var = "reads")
scaled <- cbind(temp[,1], scale(temp[,-1], center = F, scale = T))

# merge back to present
temp <- melt(scaled, id.vars = "Sample", variable.name = "OTU", value.name = "z.reads")
present <- present[temp, , on = .(Sample, OTU)][!is.na(reads),]

## Filter too rare MOs -------------------------------------------------------------------------------------
# only keep those that have more or equal then 7 observations
length(unique(levels(factor(present$ID)))) #28814
hist(present$n.obs, breaks = seq(0, max(present$n.obs) + 10, 1))

# decided based on histogram, n.obs reaches plateau at 7
present <- present[n.obs >= 7,]
length(unique(levels(factor(present$ID)))) #18471

#saveRDS(present, paste0("./Objects/MO_present_", Sys.Date(), ".rds"))

## Randomization approach --------------------------------------------------------------------------------
# For each ID, shuffle the peak intensity along the x-axis and run regression 999 times

# this is actually a mini version of the full function.
# run script 5.1_MO_randomize.R on Compute Canada

# slope.random <- ddply(x, .(ID), function(x) {
# 
#   out <- foreach(i = 1:999, .combine = rbind) %dopar% {
#   # set random iteration seed
#   set.seed(i)
#   # extract data
#   df <- data.frame(x = x$travel.time_d, y = x$z.peak.int)
#   # shuffle y over x
#   df$y <- sample(df$y, replace = F, size = nrow(df))
#   # do linear regression
#   lin <- lm(df$y ~ df$x)
#   # extract run and slope
#   out <- data.frame(run = i, slope = lin$coefficients[2])
#   return(out)
#   }
#   
#   return(out)
# }, .parallel = T)

# read in original present data frame
present <- readRDS("./Objects/MO_present_2024-02-10.rds") %>% setDT()
# read in results
out <- readRDS("./Objects/MO_randomization_output_2024-02-10.rds") %>% setDT()
ci <- readRDS("./Objects/MO_randomization_CI_2024-02-10.rds") %>% setDT()
ci <- ci[order(ID),]

# First filter, calculate the number of reshuffling that worked for a specific ID
out <- out %>% distinct() %>% setDT()
# count the number of runs that did not work
out[is.na(slope), n.na := .N, by = .(ID)]
out[is.na(n.na), n.na := 0]
# count the number of runs that did work
out[!is.na(slope), n.ok := .N, by = .(ID)]
out[is.na(n.ok), n.ok := 0]
# summarize
check <- out[, .(n.na = unique(n.na),
        n.ok = unique(n.ok)), by = .(ID)]
nrow(check) == length(unique(out$ID))

# IDs that did not work
not.working <- unique(check[n.ok == 0,]$ID)
length(not.working) # 532

# constant but above 0
const <- unique(ci[sd == 0 & mean > 0,]$ID)
length(const) # 1387

# remove those IDs that did not work
present <- present[!(ID %in% not.working),]


# Plot a few examples
set.seed(3)
exmpl <- sample(unique(present$ID), 20)

for(i in 1:length(exmpl)){
  # correlation plot
  cor <- ggplot(present[ID == exmpl[i],], aes(x = travel.time_d, y = z.reads)) +
    theme_bw() +
    geom_smooth(method = "lm") +
    geom_point()
  # histogram of p-values
  his <- ggplot(out[ID == exmpl[i],], aes(x = slope)) +
    theme_bw() +
    geom_histogram(fill = "white", colour = "gray50") +
    geom_vline(xintercept = ci[ID == exmpl[i] & vars == "slope",]$mean,
               colour = "gray20", linetype = "solid") +
    geom_vline(xintercept = ci[ID == exmpl[i] & vars == "slope",]$upper.bound,
               colour = "tomato", linetype = "dashed") +
    geom_vline(xintercept = ci[ID == exmpl[i] & vars == "slope",]$lower.bound,
               colour = "royalblue", linetype = "dashed")
  
  # histogram of R2
  r2 <- ggplot(out[ID == exmpl[i],], aes(x = r2)) +
    theme_bw() +
    geom_histogram(fill = "white", colour = "gray50")
  # histogram of p-values
  pval <- ggplot(out[ID == exmpl[i],], aes(x = p.val)) +
    theme_bw() +
    geom_histogram(fill = "white", colour = "gray50")
  
  # arranged
  arr <- ggarrange(ggarrange(cor, his, align = "v", nrow = 2), ggarrange(r2, pval, align = "hv"), nrow = 2, 
                   heights = c(1,0.5))
  ggsave(paste0("./Figures/Randomization/MO_", exmpl[i],".png"), arr,
         width = 10, height = 15, units = "cm")
}

# Apply to all MOs -------------------------------------------------------------------------------------
# for applying function see script 5.2

model.ls <- readRDS("./Objects/MO_model.ls_binned_2024-02-10.rds")
model.df <- readRDS("./Objects/MO_model.df_binned_2024-02-10.rds") %>% setDT()
pks.df <- readRDS("./Objects/MO_pks.df_binned_2024-02-10.rds") %>% setDT()

# Clean --------------------------------------------------------------------------------------------
model.df[is.na(p.val), AU := "n.s."]

nrow(model.df[is.na(p.val),]) * 100 / length(model.ls) # we loose 20.1%

model.df[model == "y ~ x", model.type := "LM"]
model.df[model == "y ~ poly(x, 3)" | model == "y ~ poly(x, 2)", model.type := "POLY"]
model.df[model == "y ~ s(x)", model.type := "GAM"]
model.df[is.na(model.type),]

# count numbers of model by type
model.df[, .(n = .N), by = .(model.type)] # 2237 NA models

model.df <- model.df %>% separate(col = "ID", into = c("year","season","dna.type","OTU"), sep = "[.]", remove = F) %>% setDT()

# Add randomization --------------------------------------------------------------------------------
model.df[ci[vars == "slope",], c("rand.ci.up", "rand.ci.low") := list(i.upper.bound, i.lower.bound), on = .(ID)]

length(unique(model.df$ID)) # 17939 models

# direction is defined 1 == positive, 2 == negative
# see how many pass the randomization filter
model.df[direction == 1 & (lm.slope > rand.ci.up), r.test := TRUE]
model.df[direction == 2 & (lm.slope < rand.ci.low), r.test := TRUE]
model.df[is.na(r.test), r.test := FALSE]

model.df[, .(n = .N), by = .(r.test, dna.type)] # 11308 DNA and 3968 RNA pass the test

# see how many pass p-value filter
model.df[p.val < 0.05, .(n = .N)] # 1949
model.df[p.val < 0.01, .(n = .N)] # 672

# see how many pass both
model.df[p.val < 0.05 & r.test == TRUE, .(n = .N), by = .(dna.type)] # D 1483, R 451
model.df[p.val < 0.01 & r.test == TRUE, .(n = .N), by = .(dna.type)] # D 535, R 131

# assign flags
model.df[, loose.filt := ifelse(r.test == TRUE, TRUE, FALSE)]
model.df[, cons.filt := ifelse(r.test == TRUE & p.val < 0.05, TRUE, FALSE)]

# Define a few truly flat models ------------------------------------------------------------------------
model.df[cons.filt == TRUE & sd.y == 0, ] # none
model.df[sd.y == 0, AU := "flat"]

# classify not significant
model.df[is.na(p.val), ns.s := "n.s."]
#model.df[p.cor == "remove", ns.s := "n.s."]

## Add ranges -----------------------------------------------------------------------------------------------
# Molecules that occur along the whole continuum
limits <- present[reads != 0, .(min.x = min(travel.time_d, na.rm = T),
            max.x = max(travel.time_d, na.rm = T)), by = .(ID)]
limits <- limits[ID %in% names(model.ls),]

model.df <- model.df[limits, , on = .(ID)] #c(".id==ID")

# model.df[min.x <= -15 & max.x >= 5, range := "entire"]
# # Molecules that only occur in freshwater
# model.df[min.x >= -15, range := "freshwater"]
# # Molecules that do not occur in estuary
# model.df[min.x <= -15 & max.x < 6, range := "soil-fresh"]


## Define model.type ------------------------------------------------------------------------
model.df[is.na(model.type), model.type := "none"]
# all defined?
model.df[is.na(model.type),] # yes

## Classify models ----------------------------------------------------------------------------------------
# Classify LMs
# classification based on their slope direction
model.df[is.na(AU) & model.type == "LM" & direction == 1, AU := "increase"] #& p.val < 0.05
model.df[is.na(AU) & model.type == "LM" & direction == 2, AU := "decrease"] #& p.val < 0.05
model.df[is.na(AU) & model.type == "LM" & lm.slope == 0, AU := "flat"]

#sanity check
model.df[model.type == "LM" & is.na(AU),]

# One peakers -------------------------------------------------------------------------------------------
# Add a few variables for classification of patterns
model.df[, num.range := max.x - min.x]
model.df[, range.cntr := (min.x + max.x) / 2]
model.df[, buffer.cntr := num.range / 6]
model.df[, buffer.min := range.cntr - buffer.cntr]
model.df[, buffer.max := range.cntr + buffer.cntr]

# Merge model vars to peaks data frame
df <- pks.df[model.df, , on = .(ID)]
df <- df[is.na(AU), ]

# Apply classification settings
df[is.na(AU) & no.peak == 1 &
     peak.x >= buffer.min & peak.x <= buffer.max, AU := "unimodal"]
df[is.na(AU) & no.peak == 1 & peak.x < range.cntr & model.type != "LM", AU := "non-linear decrease"]
df[is.na(AU) & no.peak == 1 & peak.x > range.cntr & model.type != "LM", AU := "non-linear increase"]
# sanity check, one peakers are classified?
df[is.na(AU) & no.peak == 1,]
# which were identified?
unique(df$AU)

# merge into model.df
model.df[ID %in% df[AU == "unimodal",]$ID, AU := "unimodal"]
model.df[ID %in% df[AU == "non-linear decrease",]$ID, AU := "non-linear decrease"]
model.df[ID %in% df[AU == "non-linear increase",]$ID, AU := "non-linear increase"]


# Check how many for each AU
model.df[, .(n = .N), by = .(AU)]

# Bi-peakers ---------------------------------------------------------------------------------------------------
df <- pks.df[model.df, , on = .(ID)]
df <- df[is.na(AU), ]

# Classify two-peakers ------------------------------------------------------------------------------------------
## Quadratic ----------------------------------------------------------------------------------------------------
temp <- df[is.na(AU) & no.peak == 2 & peak == 1 & peak.x <= (min.x + buffer.cntr), ]$ID
# 1st peak is at the beginning of x-axis
temp <- df[(ID %in% temp) & peak == 2 & peak.x >= (max.x - buffer.cntr),]$ID
# 2nd peak is at the end of x-axis
pks.diff <- dcast(df[ID %in% temp,], ID ~ peak, value.var = "peak.y", subset = .(closest == "above"))
pks.diff[, diff := `1`-`2`]
# when the difference is < 0, then the first peak is smaller than second = increasing
# when the difference is > 0, then the first peak is bigger than second = decreasing
model.df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing quadratic"]
model.df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing quadratic"]
df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing quadratic"]
df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing quadratic"]

# Bimodal ------------------------------------------------------------------------------------------------------
pks.diff <- dcast(df[is.na(AU) & no.peak == 2, ], ID ~ peak, value.var = "peak.y", subset = .(closest == "above"))
pks.diff[, diff := `1` - `2`]

model.df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing bimodal"]
model.df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing bimodal"]
df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing bimodal"]
df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing bimodal"]

# sanity check, all two peakers are classified?
df[is.na(AU) & no.peak == 2,]

# how many peaks max?
max(pks.df$no.peak) # 5

# what models were selected?
unique(model.df$AU)

# Classify multi-peakers ------------------------------------------------------------------------------------------
pks.diff <- dcast(df[is.na(AU), ], ID ~ peak, value.var = "peak.y", subset = .(closest == "above"))
pks.diff <- melt(pks.diff, id.vars = "ID", 
                 measure.vars = c("1","2","3","4","5"), 
                 variable.name = "peak", value.name = "y" ) %>% arrange(ID)
pks.diff <- pks.diff[!is.na(y),]

# calcualte the difference between the first and last peak
opposite.diff <- ddply(pks.diff, .(ID), function(x){
  data.frame(max.diff = x[x$peak == 1,]$y - x[which.max(x$peak),]$y)
}) %>% setDT()

opposite.diff[max.diff < 0, AU := "non-linear increase"]
opposite.diff[max.diff > 0, AU := "non-linear decrease"]
opposite.diff[is.na(AU),] # none

model.df[ID %in% unique(opposite.diff[AU == "non-linear increase",]$ID), AU := "non-linear increase"]
model.df[ID %in% unique(opposite.diff[AU == "non-linear decrease",]$ID), AU := "non-linear decrease"]

#sanity check
model.df[is.na(AU),] # none, all classified

# add info whether non-linear models have humps
cast.pks <- dcast(pks.df[no.peak != 0,], ID + peak ~ closest, value.var = "infp.slope")
model.df[ID %in% cast.pks[!is.na(above) & !is.na(below),]$ID, hump := TRUE]
model.df[is.na(hump), hump := FALSE]

# identify those that are significant
model.df[!is.na(ns.s), c.ns.s := "n.s."]
model.df[is.na(ns.s) & loose.filt == F, ns.s := "n.s."]
model.df[is.na(c.ns.s) & cons.filt == F, c.ns.s := "n.s."]
model.df[is.na(ns.s) & loose.filt == T, ns.s := "sig."]
model.df[is.na(c.ns.s) & cons.filt == T, c.ns.s := "sig."]
model.df[AU == "no.p" | AU == "stable", ns.s := "n.s."]
model.df[AU == "stable" | AU == "flat",]
# sanity check
model.df[is.na(ns.s),] #none
model.df[is.na(c.ns.s),] #none

model.df[model.type == "none", AU := "flat"]

# Save results ----------------------------------------------------------------------------------------------------
unique(model.df$AU)
model.df[, AU := factor(AU, levels = c("n.s.", "flat",
                                       "increase","non-linear increase", "decreasing quadratic","decreasing bimodal",
                                       "unimodal",
                                       'decrease','non-linear decrease', 'increasing quadratic','increasing bimodal'))]

model.df <- model.df %>% dplyr::select(ID:p.val, min.x:buffer.max, rand.ci.up:cons.filt, model.type, hump, AU, ns.s, c.ns.s)

write.table(model.df,paste0("./Data/Summary/MO_patterns.classified_", Sys.Date(), ".csv"), sep = ",", row.names = F)



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
#   [1] ade4_1.7-23         ape_5.8-1           vegan_2.7-2         permute_0.9-10     
# [5] kableExtra_1.4.0    cowplot_1.2.0       phyloseq_1.54.2     BiocManager_1.30.27
# [9] plotly_4.12.0       gridExtra_2.3       doMC_1.3.8          mgcv_1.9-4         
# [13] nlme_3.1-168        rstatix_0.7.3       doParallel_1.0.17   iterators_1.0.14   
# [17] foreach_1.5.2       ggpubr_0.6.3        plyr_1.8.9          lubridate_1.9.5    
# [21] forcats_1.0.1       stringr_1.6.0       dplyr_1.2.0         purrr_1.2.1        
# [25] readr_2.2.0         tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.2      
# [29] tidyverse_2.0.0     data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] rlang_1.1.7         magrittr_2.0.4      otel_0.2.0          compiler_4.5.2     
# [5] systemfonts_1.3.1   vctrs_0.7.1         reshape2_1.4.5      pkgconfig_2.0.3    
# [9] crayon_1.5.3        fastmap_1.2.0       backports_1.5.0     XVector_0.50.0     
# [13] labeling_0.4.3      CoprManager_0.5.8   rmarkdown_2.30      tzdb_0.5.0         
# [17] ragg_1.5.0          xfun_0.56           jsonlite_2.0.0      biomformat_1.38.0  
# [21] rhdf5filters_1.22.0 Rhdf5lib_1.32.0     broom_1.0.12        cluster_2.1.8.2    
# [25] R6_2.6.1            stringi_1.8.7       RColorBrewer_1.1-3  car_3.1-5          
# [29] Rcpp_1.1.1          Seqinfo_1.0.0       knitr_1.51          IRanges_2.44.0     
# [33] Matrix_1.7-4        splines_4.5.2       igraph_2.2.2        timechange_0.4.0   
# [37] tidyselect_1.2.1    rstudioapi_0.18.0   abind_1.4-8         codetools_0.2-20   
# [41] lattice_0.22-9      Biobase_2.70.0      withr_3.0.2         S7_0.2.1           
# [45] evaluate_1.0.5      survival_3.8-6      xml2_1.5.2          Biostrings_2.78.0  
# [49] pillar_1.11.1       carData_3.0-6       stats4_4.5.2        generics_0.1.4     
# [53] S4Vectors_0.48.0    hms_1.1.4           scales_1.4.0        glue_1.8.0         
# [57] lazyeval_0.2.2      tools_4.5.2         ggsignif_0.6.4      rhdf5_2.54.1       
# [61] grid_4.5.2          Formula_1.2-5       cli_3.6.5           textshaping_1.0.4  
# [65] viridisLite_0.4.3   svglite_2.2.2       gtable_0.3.6        digest_0.6.39      
# [69] BiocGenerics_0.56.0 htmlwidgets_1.6.4   farver_2.1.2        htmltools_0.5.9    
# [73] multtest_2.66.0     lifecycle_1.0.5     httr_1.4.8          MASS_7.3-65 