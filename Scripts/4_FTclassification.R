## -------------------------------------------------------------------------
##
## Script name: 4_FTclassification.R
##
## Purpose of script: This script is the main script for FT analyses.
##                    We identify reactive versus unreactive molecular formulae
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
##         4.1 and 4.2 and then, their outputs are continued to be
##         analyzed in this script.
##         The script is written in this way, because some sections
##         were run on a high-performance computer provided by
##         Compute Canada.
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
  "data.table", "stringr", "tidyverse", # wrangling
  "plyr",
  "rstatix", "mgcv", #stats 
  "doMC", "foreach", # parallel
  'gridExtra', 'ggpubr','plotly') # plotting

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

cross <- read.csv("./Data/FT/crosstable2015_cor.csv", sep =",", stringsAsFactors = F)
cross16 <- read.csv("./Data/FT/crosstable2016_cor.csv", sep =",", stringsAsFactors = F)

# # read in matched mother files
# # we will not be publishing all variables we have sampled
# # make a selection
# meta15 <- read.csv("./Data/Summary/2015_metadata_upflacc.csv", sep = ",", stringsAsFactors = F) %>% setDT()
# meta16 <- read.csv("./Data/Summary/2016_metadata_upflacc.csv", sep = ",", stringsAsFactors = F)%>% setDT()
# 
# meta15.all <- read.csv("./Data/Summary/2015_metadata.csv", sep = ",", stringsAsFactors = F) %>% setDT()
# meta16.all <- read.csv("./Data/Summary/2016_metadata.csv", sep = ",", stringsAsFactors = F) %>% setDT()
# 
# # 2015
# meta15 <- meta15[meta15$sample.name %in% meta15.all$sample.name]
# 
# meta15 <- left_join(meta15.all, meta15 %>% dplyr::select(samples, elev_m:snapped_long_tm5), by = "samples")
# meta15 <- meta15 %>% dplyr::select(samples:sample.type, sample.type.year, Season,
#                          lat:long, temp.water.ysi,do.cor.perc.ysi, ph.ysi, DOC, bact.abundance, elev_m:snapped_long_tm5)
# write.table(meta15, "./Data/Summary/2015_metadata_selection.csv", sep = ",", row.names = F)
# 
# # 2016
# meta16 <- meta16[meta16$sample.name %in% meta16.all$sample.name]
# 
# meta16 <- left_join(meta16.all, meta16 %>% dplyr::select(samples, elev_m:snapped_long_tm5), by = "samples")
# meta16 <- meta16 %>% dplyr::select(samples:sample.type, sample.type.year, Season,
#                                    lat:long, temp.water.ysi,do.cor.perc.ysi, ph.ysi, DOC, bact.abundance, elev_m:snapped_long_tm5)
# write.table(meta16, "./Data/Summary/2016_metadata_selection.csv", sep = ",", row.names = F)

# read in mother files
meta15 <- read.csv("./Data/Summary/2015_metadata_selection.csv", sep = ",", stringsAsFactors = F)
meta16 <- read.csv("./Data/Summary/2016_metadata_selection.csv", sep = ",", stringsAsFactors = F)

meta.all <- read.csv("./Data/Summary/meta_file_traveltime_2024-02-10.csv", sep = ",", stringsAsFactors = F) %>% setDT()


# Format data ------------------------------------------------------------
# extract only a few samples that we have FT for
meta <- meta.all[sample.name %in% c(meta15$sample.name, meta16$sample.name),] %>% unique()

# join
meta <- left_join(meta,
                  rbind(meta15, meta16) %>% dplyr::select(sample.name, colnames(meta15)[!(colnames(meta15) %in% colnames(meta.all))]),
                  by = "sample.name")

meta[, (n = .N), by = .(year)]

# calculate some more FT variables
cross.tabs <- list(cross, cross16)

# Add FT variables -----------------------------------------------------------------------------------------

ios <- function(x){
  setDT(x)
  x[, molweight := (12.0107*C) + (1.00784*H) + (14.0067*N) + (15.999*O) + (32.065*S)]
  x[, IOS := F]
  x[HC <= 1.3 & HC >= 1.04 & OC <= 0.62 & OC >= 0.42 & molweight <= 388 & molweight >= 332, IOS := T]
  x[HC <= 1.3 & HC >= 1.04 & OC <= 0.62 & OC >= 0.42 & molweight <= 548 & molweight >= 446, IOS := T]
  return(x)
}
cross.tabs <- lapply(cross.tabs, ios)


# Merge and clean --------------------------------------------------------------------------------------
# melt normal samples into long format
cross.tabs[[1]] <- melt.data.table(cross.tabs[[1]], id.vars = c(1:12,15:18,21:26,122:123),
                                   measure.vars = 33:121,
                                   variable.name = "sample.name", value.name = "peak.int")
cross.tabs[[2]] <- melt.data.table(cross.tabs[[2]], id.vars = c(1:12,15:18,21:26,126:127),
                                   measure.vars = 33:125,
                                   variable.name = "sample.name", value.name = "peak.int")

# # merge together
all.df <- bind_rows(cross.tabs)

## Sample renaming ---------------------------------------------------------------------------------------
# some cleaning necessary
# change sample names from . to -
all.df[str_detect(sample.name, "^RO"), sample.name := str_replace(sample.name, "[.]", "-")]
meta[str_detect(sample.name, "^PR"), sample.name := str_replace(sample.name, "[-]", ".")]

all.df[str_detect(sample.name, "^L32"), sample.name := str_replace(sample.name, "[.]", "_")]

# RO1-10
meta[sample.name == "RO1-10D", sample.name := "RO1-10"]
all.df[sample.name == "RO1-10P", sample.name := "RO1-10"]

# RO2-22.r2
meta[sample.name == "RO2-22", sample.name := c("RO2-22","RO2-22.r2")]

# omit estuaries
est <- meta[sample.type.year == "Estuary",]$sample.name
all.df <- all.df[!(sample.name %in% est),]
meta <- meta[sample.type.year != "Estuary",]

#write.table(meta, paste0("./Data/Summary/meta_FTdataset_", Sys.Date(), ".csv"), sep = ",", dec = ".", row.names = F)

# ggplot(meta[sample.type.year != "Lake" & sample.type.year != "Riverine Lakes",], aes(x = sample.type.year, y =wrt_min)) +
#   geom_boxplot()

# merge with meta
all.df <- all.df[meta, , on = .(sample.name)] 

setDT(all.df)

## Define habitats -----------------------------------------------------------------------------------------
all.df$sample.type.year <- factor(all.df$sample.type.year, levels = c(#"Soil","Sediment",
  "Groundwater","Soilwater",#"Hyporheic", "Wellwater",
  "Stream", "Tributary",
  #"IslandLake",
  "Upriver", "Reservoirs",#"Deep",
  "Downriver",
 #"Estuary",
  "Headwater Ponds", "Riverine Lakes", "Lake"))

# save coordinates
out <- all.df %>% dplyr::select(sample.name, Season, year, sample.type.year, lat, long)
out <- out %>% distinct()
#write.table(out, "./Output/coordinates_FT.csv", sep = ",", row.names = F)

# df <- all.df %>% dplyr::select(sample.name, sample.type.year,
#                         DOC, lability, bact.abundance, TP, TN, temp.water.ysi, do.cor.perc.ysi) %>% distinct()
# 
# df <- melt(df, id.vars = c("sample.name","sample.type.year"))
# 
# ggplot(df) +
#   theme_bw()+
#   geom_boxplot(aes(y = value, x = sample.type.year)) +
#   facet_wrap(.~variable, scale = "free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# setDT(df)
# options(scipen =9999)
# df[, .(mean = mean(value, na.rm = T)), by = .(variable,  sample.type.year)]


# Table S2 ----------------------------------------------------------------
# Number of samples used to create spatial models of microbial OTUs and DOM MF per campaign
# make sure we only count those we are using
temp <- meta[meta$sample.name %in% all.df$sample.name,] %>% setDT()
# calculate sample size by campaign
temp <- temp[,.(n = .N), by = .(year, Season)] %>% arrange(year, Season)

write.table(temp,"./Data/Summary/FT_samplesize.csv", sep = ",", row.names = F)

## Remove outliers -------------------------------------------------------------------------------------------
# During plotting (later step in this script), a few outliers were detected
# We apply a function to remove an outlier observation per group
# We identify a group by: year and molecular formula
# If a observation was identified as an outlier, it will be removed
# However, if there are several points of outliers they will be kept
all.df[, ID := paste(year, molecular.formula, sep = "_")]

# detect extreme group outliers, except zeros
all.df[peak.int > 0, ext_outlier := is_extreme(peak.int), by = .(Year, Season, molecular.formula)]

# overwrite all zeros with FALSE
all.df[peak.int == 0, ext_outlier := FALSE]
# check how many outliers
nrow(all.df[ext_outlier == TRUE,]) * 100 / nrow(all.df) # 1.05% outliers
# remove extreme outliers
clean.df <- all.df[ext_outlier == TRUE, peak.int := 0]

## Separate absent and present MF ---------------------------------------------------------------------------
# Define absence
clean.df[, ID := paste(Year, Season, molecular.formula, sep = "_")]
clean.df[, n := .N, by = .(ID)] # number of samples per category
clean.df[, n.obs := nrow(.SD[peak.int > 0,]), by = .(ID)] # number of samples with actual observation

#initiate absence column
clean.df[n.obs == 0, PA := "Absent", by = .(ID)]

levels(factor(clean.df$sample.type.year))

#how many molecular formulae that are absent?
length(levels(factor(clean.df[PA == "Absent",]$molecular.formula))) # 2283
length(levels(factor(clean.df$molecular.formula))) # 13876
# percentage lost
length(levels(factor(clean.df[PA == "Absent",]$molecular.formula))) * 100 / length(levels(factor(clean.df$molecular.formula)))
# 16.5% absent

# Overwrite the rest with present
clean.df[is.na(PA), PA := "Present"]

levels(factor(clean.df[PA == "Present",]$sample.type.year))
levels(factor(clean.df[PA == "Absent",]$sample.type.year))
absent <- clean.df[PA == "Absent",]

# extract only df with present molecular formulae
present <- clean.df[PA == "Present",]
levels(factor(present$sample.type.year))

## Define travel time for non estimated habitats ----------------------------------------------------------
present[, travel.time_d := travel.time_days]
present[sample.type.year == "Soilwater", travel.time_d := -50]
present[sample.type.year == "Groundwater", travel.time_d := -50]

## Transform variables and remove zero observations --------------------------------------------------------
present[peak.int == 0, peak.int := NA]

# look at an example distribution along travel time
ggplot(present[molecular.formula == "C10H10O10N0S0",], 
       aes(x =travel.time_d, y = peak.int, colour = sample.type.year)) +
  geom_point()

# arrange so that the data frame follows travel time
# ID is year_campaign_molecular.formula
present <- present %>% group_by(ID) %>% dplyr::arrange(travel.time_d) %>% setDT()

## Z-Scale -------------------------------------------------------------------------------------------------
temp <- dcast(present, sample.name ~ molecular.formula, value.var = "peak.int")
# randomize
#temp <- temp[, (colnames(temp)[-1]) := lapply(.SD, replace_na, 0), .SDcols = colnames(temp)[-1]]
scaled <- cbind(temp[,1], scale(temp[,-1], center = F, scale = T))

# merge back to present
temp <- melt(scaled, id.vars = "sample.name", variable.name = "molecular.formula", value.name = "z.peak.int")
present <- present[temp, , on = .(sample.name, molecular.formula)][!is.na(peak.int),]

# remove Autumn -----------------------------------------------------------------------------------------
present <- present[Season != "Autumn",]

## Filter too rare MF -------------------------------------------------------------------------------------
# Only keep those with more than 8 observations
length(unique(levels(factor(present$ID)))) #42484
hist(present$n.obs, breaks = seq(0, max(present$n.obs) + 10, 1))
# histogram is hard, three major peaks
present <- present[n.obs >= 8,]
length(unique(levels(factor(present$ID)))) #36193

#saveRDS(present, paste0("./Objects/FT_present_", Sys.Date(), ".rds"))

## Randomization approach --------------------------------------------------------------------------------
# For each ID, shuffle the peak intensity along the x-axis and run regression 999 times

# this is actually a mini example of the full function.
# run script 4.1_FT_randomize.R on Compute Canada

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
present <- readRDS("./Objects/FT_present_2024-02-10.rds") %>% setDT()
# read in results
out <- readRDS("./Objects/FT_randomization_output_2024-02-10.rds") %>% setDT()
ci <- readRDS("./Objects/FT_randomization_CI_2024-02-10.rds") %>% setDT()
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
length(not.working) # 0

# constant but above 0
const <- unique(ci[sd == 0 & mean > 0,]$ID)
length(const) # 0

# Plot a few examples
set.seed(3)
exmpl <- sample(unique(present$ID), 20)

for(i in 1:length(exmpl)){
  # correlation plot
  cor <- ggplot(present[ID == exmpl[i],], aes(x = travel.time_d, y = z.peak.int)) +
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
  ggsave(paste0("./Figures/Randomization/", exmpl[i],".png"), arr,
         width = 10, height = 15, units = "cm")
}

# For model application refer to script 4.2
# Classification ------------------------------------------------------------------------------------
# read objects
model.ls <- readRDS("./Objects/FT_model.ls_binned_2024-02-11.rds")
model.df <- readRDS("./Objects/FT_model.df_binned_2024-02-11.rds") %>% setDT()
pks.df <- readRDS("./Objects/FT_pks.df_binned_2024-02-11.rds") %>% setDT()

# Plot examples ------------------------------------------------------------------------------------
# remove NA models (not enough observations = none)

# # plot tests
# plots <- lapply(model.ls[1:10], function(x){
#   newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
#   newd$y.pred <- predict(x$model, newd)
#   p <- ggplot()+
#     theme_bw() +
#     geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
#     geom_point(aes(x = x$binned$x, y = x$binned$y,
#                    fill = factor(x$binned$n, levels = seq(min(x$binned$n),
#                                                             max(x$binned$n), 1))),
#                size = 5, shape = 21) +
#     scale_fill_brewer(palette = "RdBu", direction = -1, name = "Sample size") +
#     geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#     labs(title = paste("Model = ", x$model.df$model, "\nDirection = ", x$model.df$direction),
#          x = "Travel time (d)",
#          y = "Peak intensity")
# 
#   if(!is.null(x$pks.df)){
#     if(x$pks.df$no.peak[1] > 0){
#     p <- p +
#       geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", size = 4, shape = 21) +
#       geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2) +
#       scale_colour_manual(values = c("tomato","skyblue"), name = "Inflection point")
#   }}
#   return(p)
# })
# 
# # save plots
# for(i in 1:10){
#   ggsave(paste0("./Figures/FT_lm_fit/", names(plots)[i],".png"), plots[[i]])
# }

# Extract raw data
# raw.df <- ldply(model.ls, function(x){
#   x$raw
# }) %>% setDT()

# Clean --------------------------------------------------------------------------------------------
model.df[is.na(p.val), AU := "n.s."]

nrow(model.df[is.na(p.val),]) * 100 / length(model.ls) # we loose 0.58%

model.df[model == "y ~ x", model.type := "LM"]
model.df[model == "y ~ poly(x, 3)" | model == "y ~ poly(x, 2)", model.type := "POLY"]
model.df[model == "y ~ s(x)", model.type := "GAM"]
model.df[is.na(model.type),]

# count numbers of model by type
model.df[, .(n = .N), by = .(model.type)]

model.df <- model.df %>% separate(col = "ID", into = c("year","season","MF"), sep = "[_]", remove = F) %>% setDT()

# Add randomization --------------------------------------------------------------------------------
model.df[ci[vars == "slope",], c("rand.ci.up", "rand.ci.low") := list(i.upper.bound, i.lower.bound), on = .(ID)]

length(unique(model.df$ID)) # 36193 models

# direction is defined 1 == positive, 2 == negative
# see how many pass the randomization filter
model.df[direction == 1 & (lm.slope > rand.ci.up), r.test := TRUE]
model.df[direction == 2 & (lm.slope < rand.ci.low), r.test := TRUE]
model.df[is.na(r.test), r.test := FALSE]

model.df[, .(n = .N), by = .(r.test)] # 35312 pass

# see how many pass p-value filter
model.df[p.val < 0.05, .(n = .N)] # 16073
model.df[p.val < 0.01, .(n = .N)] # 9477

# see how many pass both
model.df[p.val < 0.05 & r.test == TRUE, .(n = .N)] # 15905
model.df[p.val < 0.01 & r.test == TRUE, .(n = .N)] # 9384

# assign flags
model.df[, loose.filt := ifelse(r.test == TRUE, TRUE, FALSE)]
model.df[, cons.filt := ifelse(r.test == TRUE & p.val < 0.05, TRUE, FALSE)]

# Define a few truly flat models ------------------------------------------------------------------------
model.df[cons.filt == TRUE & sd.y == 0, ] # none
model.df[sd.y == 0, AU := "flat"]

# Classification ------------------------------------------------------------------------------------
# classify not significant
model.df[is.na(AU) & is.na(p.val), ns.s := "n.s."]

## Add ranges -----------------------------------------------------------------------------------------------
# Molecules that occur along the whole continuum
limits <- present[peak.int != 0, .(min.x = min(travel.time_d, na.rm = T),
                                max.x = max(travel.time_d, na.rm = T)), by = .(ID)]
limits <- limits[ID %in% names(model.ls),]

model.df <- model.df[limits, , on = .(ID)] #c(".id==ID")

# model.df[min.x <= -15 & max.x >= 5, range := "entire"]
# # Molecules that only occur in freshwater
# model.df[min.x >= -15, range := "freshwater"]
# # Molecules that do not occur in estuary
# model.df[min.x <= -15 & max.x < 6, range := "soil-fresh"]

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
                 measure.vars = c("1","2","3","4"), 
                 variable.name = "peak", value.name = "y" ) %>% arrange(ID)
pks.diff <- pks.diff[!is.na(y),]

# calculate the difference between the first and last peak
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
# quality filter
model.df[is.nan(p.val), c("ns.s", "c.ns.s") := list("n.s.", "n.s.")]
# sanity check
model.df[is.na(ns.s),] #none
model.df[is.na(c.ns.s),] #none

model.df[model.type == "none", AU := "flat"]


# Save classification -----------------------------------------------------------------------
model.df[is.na(p.val), AU := "n.s."]
unique(model.df$AU)
model.df[, AU := factor(AU, levels = c("n.s.", "flat",
                                       "increase","non-linear increase", "decreasing quadratic","decreasing bimodal",
                                       "unimodal",
                                       'decrease','non-linear decrease', 'increasing quadratic','increasing bimodal'))]

model.df <- model.df %>% dplyr::select(ID:p.val, min.x:buffer.max, rand.ci.up:cons.filt, model.type, hump, AU, ns.s, c.ns.s)

write.table(model.df,paste0("./Data/Summary/FT_patterns.classified_", Sys.Date(), ".csv"), sep = ",", row.names = F)





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
#   [1] plotly_4.12.0       gridExtra_2.3       doMC_1.3.8          mgcv_1.9-4          nlme_3.1-168       
# [6] rstatix_0.7.3       doParallel_1.0.17   iterators_1.0.14    foreach_1.5.2       ggpubr_0.6.3       
# [11] plyr_1.8.9          lubridate_1.9.5     forcats_1.0.1       stringr_1.6.0       dplyr_1.2.0        
# [16] purrr_1.2.1         readr_2.2.0         tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.2      
# [21] tidyverse_2.0.0     data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6       htmlwidgets_1.6.4  lattice_0.22-9     tzdb_0.5.0         CoprManager_0.5.8 
# [6] vctrs_0.7.1        tools_4.5.2        generics_0.1.4     pkgconfig_2.0.3    Matrix_1.7-4      
# [11] RColorBrewer_1.1-3 S7_0.2.1           lifecycle_1.0.5    compiler_4.5.2     farver_2.1.2      
# [16] codetools_0.2-20   carData_3.0-6      htmltools_0.5.9    lazyeval_0.2.2     Formula_1.2-5     
# [21] pillar_1.11.1      car_3.1-5          abind_1.4-8        tidyselect_1.2.1   digest_0.6.39     
# [26] stringi_1.8.7      labeling_0.4.3     splines_4.5.2      fastmap_1.2.0      grid_4.5.2        
# [31] cli_3.6.5          magrittr_2.0.4     broom_1.0.12       withr_3.0.2        scales_1.4.0      
# [36] backports_1.5.0    timechange_0.4.0   httr_1.4.8         ggsignif_0.6.4     hms_1.1.4         
# [41] viridisLite_0.4.3  rlang_1.1.7        Rcpp_1.1.1         glue_1.8.0         jsonlite_2.0.0    
# [46] rstudioapi_0.18.0  R6_2.6.1