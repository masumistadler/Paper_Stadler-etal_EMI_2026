## -------------------------------------------------------------------------
##
## Script name: 5.1_MO_randomize.R
##
## Purpose of script: This is a randomization process done to evaluate
##                    which model's slopes are higher than purely by chance.
##
## Author: Masumi Stadler
##
## Date Finalized: 2024-02-10
##
## Copyright (c) Masumi Stadler, 2026
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes:      This script is part of the 5_MOclassification.R script,
##                    however, a dedicated separate script was written
##                    to run the code on a high performance computer.
##
## -------------------------------------------------------------------------

## Use R project with regular scripts, all paths are relative 

# Server set-up -----------------------------------------------------------
## Working directory is set from where the job is submitted
## Load library path, if on a server
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )

# R-setup -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
pckgs <- list("data.table", "tidyverse", # wrangling
              "plyr",
              "doMC","foreach") # change as needed

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
# Server version
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
registerDoMC(cores = cores)

## Personal version
# detectCores(); registerDoMC(); getDoParWorkers()
# numCores <- detectCores()
# cl <- makeCluster(numCores, type = "FORK")

# Read in data -----------------------------------------------------------
present <- readRDS("./Objects/MO_present_2024-02-10.rds") %>% setDT()

# Run function -----------------------------------------------------------

random <- ddply(present, .(ID), function(x) {
  
  out <- foreach(i = 1:999, .combine = rbind) %dopar% {
    # set random iteration seed
    set.seed(i)
    # extract data
    df <- data.frame(x = x$travel.time_d, y = x$z.reads)
    # shuffle y over x, needs to be with replacement for bootstrapping
    df$y <- sample(df$y, replace = T, size = nrow(df))
    # do linear regression
    lin <- try(lm(df$y ~ df$x), silent = T)
    # extract run and slope
    if (!inherits(lin, "try-error") & !is.na(lin$coefficients[[2]])) {
    out <- data.frame(run = i, slope = lin$coefficients[2], r2 = summary(lin)$r.squared,
                      p.val = summary(lin)$coefficients[2,4])
    } else {
    out <- data.frame(run = 1, slope = NA, r2 = NA,
                      p.val = NA)
    }
    return(out)
  }
  
  return(out)
}, .parallel = T)

# Save intermediate output -----------------------------------------------------------------------------
saveRDS(random,paste0("./Objects/MO_randomization_output_", Sys.Date(), ".rds"))
#write.table(random,"./Objects/MO_randomization_output.csv", sep = ',', dec = ".", row.names = F)

# Calculate confidence interval ------------------------------------------------------------------------
setDT(random)
random <- melt(random, id.vars = 1:2, measure.vars = 3:5, variable.name = "vars", value.name = "value")

df <- random[, .(mean = mean(value, na.rm =T),
                 sd = sd(value, na.rm = T),
                 n = .N), by = .(ID, vars)]

df[, c("se","alpha", "df") := list(sd / sqrt(n),
                                   0.05,
                                   n - 1)]
df[, t.score := qt(p=alpha/2, df = df, lower.tail = F)]
df[, margin.error := t.score * se]
df[, c("upper.bound", "lower.bound") := list(mean + margin.error,
                                             mean - margin.error)]

# Save -----------------------------------------------------------------------------------------------
saveRDS(df,paste0("./Objects/MO_randomization_CI_", Sys.Date(), ".rds"))
#write.table(df,"./Objects/MO_randomization_CI.csv", sep = ',', dec = ".", row.names = F)
