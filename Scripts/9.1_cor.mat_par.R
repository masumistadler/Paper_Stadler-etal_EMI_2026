## -------------------------------------------------------------------------
##
## Script name: 9.1_cor.mat_par.R
##
## Purpose of script: This script was written to run on a HPC due to it's
##                    computationally heavy nature.
##                    It takes the binned DOM MF and microbial OTU data and
##                    correlates them to each other within each year and
##                    season.
##
## Author: Masumi Stadler
##
## Date Finalized: 2024-03-04
##
## Copyright (c) Masumi Stadler, 2026
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes: All intermediate files are very big, hence, we are only
##        providing in the initial input data and the final cleaned
##        correlation matrix.
##
## -------------------------------------------------------------------------

## Use R project with regular scripts, all paths are relative 

# Server set-up -----------------------------------------------------------
# Working directory is set from where the job is submitted
# Load library path, if on a server
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )

# R-setup -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
pckgs <- list("plyr", "tidyverse", "data.table", # wrangling
              "doMC", "doParallel", "foreach", "doMC", "doSNOW") # statistics

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

## Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# # Load custom functions --------------------------------------------------
# funs <- list.files("./Functions", full.names = T)
# invisible(lapply(funs, source))
# 
# # Other set-up -----------------------------------------------------------
# options(scipen = 6, digits = 4) # view outputs in non-scientific notation

# Parallel environment ---------------------------------------------------
# Server version
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
registerDoMC(cores = cores)

# # Personal version
# detectCores(); registerDoMC(); getDoParWorkers()
# numCores <- detectCores()
# cl <- makeCluster(numCores, type = "FORK")
 
# Read in data -----------------------------------------------------------
# get binned values from each data set
ft.model.ls <- readRDS("./Objects/FT_model.ls_binned_2024-02-11.rds")

ft.bin <- ldply(ft.model.ls, function(x){
  x$binned[,c("x", "y")]
})

# split ID
ft.bin <- ft.bin %>% separate(col = ID, into = c("Year","Season","MF"))
setDT(ft.bin)
ft.bin <- dcast(ft.bin, Year + Season + x ~ MF, value.var = "y")

# split them into year/season bins
ft.bin.ls <- dlply(ft.bin, .(Year, Season), function(x){x %>% dplyr::select(-Year, -Season)})

# Do the same for microbial data set
mo.model.ls <- readRDS("./Objects/MO_model.ls_binned_2024-02-10.rds")

mo.bin <- ldply(mo.model.ls, function(x){
  x$binned[,c("x", "y")]
})

# split ID
mo.bin <- mo.bin %>% separate(col = ID, into = c("Year","Season","dna.type","OTU"), sep = "[.]")
mo.bin <- setDT(mo.bin)
mo.bin <- mo.bin[dna.type == "DNA",]
mo.bin <- dcast(mo.bin, Year + Season + x ~ OTU, value.var = "y")

# split them into year/season bins
mo.bin.ls <- dlply(mo.bin, .(Year, Season), function(x){x %>% dplyr::select(-Year, -Season)})

# remove columns with all NAs
mo.bin.ls <- lapply(mo.bin.ls, function(x){
  x[, colSums(is.na(x)) < nrow(x)]
})

ft.bin.ls <- lapply(ft.bin.ls, function(x){
  x[, colSums(is.na(x)) < nrow(x)]
})

# clean memory
rm(ft.model.ls, mo.model.ls, mo.bin, ft.bin, cor.mat, p.mat, x, y)


# Do correlations --------------------------------------------------------------
# Set seed for session and reproducibility of permutations
set.seed(3)

# 4 correlation matrices (year x season groups)

for(m in 1:4){
  x <- mo.bin.ls[[m]]
  y <- ft.bin.ls[[m]]
  
  # start ----
  # make sure that x-axis corresponds between two matrices
  if(nrow(x) != nrow(y)){
    x <- x[which(x$x %in% y$x),]
    y <- y[which(y$x %in% x$x),]
  }
  
  # merge the two
  z <- merge(x,y, by = "x")
  
  # remove x axis from matrix
  row.names(z) <- z$x; z$x <- NULL
  
  # create structure of matrix to fill
  cor.mat <- matrix(nrow = ncol(z), ncol = ncol(z))
  colnames(cor.mat) <- colnames(z)
  row.names(cor.mat) <- colnames(z)
  p.mat <- cor.mat # duplicate to fill with p-values
  
  for(i in 1:ncol(z)){
    for(j in 1:ncol(z)){
      
      spc <- try(cor.test(z[,i], z[,j], method = "spearman", exact = FALSE, na.action = "na.omit"), silent = T)
      
      if(!inherits(spc, 'try-error')){
        cor.mat[j,i] <- spc$estimate
        cor.mat[i,j] <- spc$estimate
        
        p.mat[j,i] <- spc$p.value
        p.mat[i,j] <- spc$p.value
      } else {
        cor.mat[i,j] <- NA
        p.mat[i,j] <- NA
      }
    }
  }
  
  out.ls <- list(cor = cor.mat, p.val = p.mat)
  saveRDS(out.ls, paste0("./Objects/cor.mat_list_",names(mo.bin.ls)[m],"_", Sys.Date(),".rds"))
  rm(out.ls, x, y, cor.mat, p.mat)
}

close(pb)
stopCluster(cl)


# Read in correlation matrices -------------------------------------------------------------------------------------------------
mat.vec <- c("./Objects/cor.mat_list_2015.Spring_2024-02-15.rds",
             "./Objects/cor.mat_list_2015.Summer_2024-02-17.rds",
             "./Objects/cor.mat_list_2016.Spring_2024-02-17.rds",
             "./Objects/cor.mat_list_2016.Summer_2024-02-18.rds")

# cleaning function for later
na_to_zero <-function(DT, pDT){
  for (j in seq_len(ncol(DT))){
    set(pDT, which(is.na(DT[[j]])), j, 0) # apply same correction to p-value dataframe
    set(DT,which(is.na(DT[[j]])), j, 0)
  }
}

ns_to_zero <-function(DT, pDT){
  for (j in seq_len(ncol(DT))){
    set(DT,which(pDT[[j]] >= 0.05), j, 0) # correct correlation matrix, remove those that are not significant
    set(pDT, which(pDT[[j]] >= 0.05), j, 0) # apply the same to p-value matrix afterwards
    
  }
}

for(i in 1:4){
  dt <- readRDS(mat.vec[i])
  dim(dt[[1]]) # 10926, 10926
  
  mfs <- row.names(dt[[1]])
  # sanity check
  any(row.names(dt[[1]]) != row.names(dt[[2]])) # should be false
  
  dt <- lapply(dt, function(x) setDT(as.data.frame(x), keep.rownames = "MF"))
  
  dt[[1]][1:10, (ncol(dt[[1]])-10):ncol(dt[[1]])]
  
  # quality check
  print(paste("Any NAs?:", any(is.na(dt[[1]]))))
  
  # check if entire row/column or single incidents
  print(paste("Any OTUs with no cor?:", any(colSums(is.na(dt[[1]])) == nrow(dt[[1]]))))
  # yes, we have OTUs that have no correlations
  print(paste("Any MFs with no cor?:", any(rowSums(is.na(dt[[1]])) == ncol(dt[[1]]))))
  # no, we don't have MFs that have no correlations
  
  # We have single incidents and entire columns
  # Fix single incidents and then remove empty columns
  na_to_zero(dt[[1]], dt[[2]])
  # sanity check
  print(paste("NAs gone?:",  any(is.na(dt[[1]]))))
  # should be FALSE
  
  # remove any columns without any correlations
  mfs <- dt[[1]]$MF
  setDF(dt[[1]][, MF := NULL], rownames = mfs)
  setDF(dt[[2]][, MF := NULL], rownames = mfs)
  dt[[2]] <- dt[[2]][,-(which(colSums(dt[[1]])==0))]
  dt[[1]] <- dt[[1]][,-(which(colSums(dt[[1]])==0))]
  
  saveRDS(dt, gsub("cor.mat_list", "clean_cor.mat_list", mat.vec[i]))
}

# Melt correlation matrix ---------------------------------------------------------------------------------------
mat.vec <- c("./Objects/clean_cor.mat_list_2015.Spring_2024-02-15.rds",
             "./Objects/clean_cor.mat_list_2015.Summer_2024-02-17.rds",
             "./Objects/clean_cor.mat_list_2016.Spring_2024-02-17.rds",
             "./Objects/clean_cor.mat_list_2016.Summer_2024-02-18.rds")

for(i in 1:4){
  dt <- readRDS(mat.vec[i])
  mfs <- row.names(dt[[1]])
  
  t <- replace(dt[[1]], lower.tri(dt[[1]], TRUE), NA)
  tp <- replace(dt[[2]], lower.tri(dt[[2]], TRUE), NA)
  setDT(t, keep.rownames = "cor_y"); setDT(tp, keep.rownames = "cor_y")
  t <- melt(t, id.vars = "cor_y", measure.vars = colnames(dt[[1]])[-1],
            variable.name = "cor_x", value.name = "rho",
            na.rm = T)
  tp <- melt(tp, id.vars = "cor_y", measure.vars = colnames(dt[[2]])[-1],
             variable.name = "cor_x", value.name = "pval",
             na.rm = T)
  
  cor.df <- merge(t, tp, by = c("cor_y", "cor_x"), all = T)
  # extract only OTU ~ MF correlations
  otu.mf <- cor.df[grep("OTU", cor_x),][grep("C", cor_y),]
  colnames(otu.mf)[1:2] <- c("MF", "OTU")
  if(nrow(otu.mf) == 0L){
    otu.mf <- cor.df[grep("OTU", cor_y),][grep("C", cor_x),]
    colnames(otu.mf)[1:2] <- c("OTU", "MF")
  }
  
  saveRDS(otu.mf, gsub("clean_cor.mat_list", "cor_otumf", mat.vec[i]))
}

# Merge all correlation matrices together ---------------------------------------------------------------------
mat.vec <- c(
  "./Objects/cor_otumf_2015.Spring_2024-02-15.rds",
  "./Objects/cor_otumf_2015.Summer_2024-02-17.rds",
  "./Objects/cor_otumf_2016.Spring_2024-02-17.rds",
  "./Objects/cor_otumf_2016.Summer_2024-02-18.rds"
)

cordf <- llply(mat.vec, function(x){
  df <- readRDS(x)
  df$Year <- substr(x, 21, 24)
  df$Season <- substr(x, 26, 31)
  return(df)
})

# merge
cordf <- bind_rows(cordf)

# save final cleaned correlation matrix
saveRDS(cordf, "./Objects/all_cors_2024-03-04.rds")
