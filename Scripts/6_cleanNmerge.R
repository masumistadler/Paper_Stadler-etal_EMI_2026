## -------------------------------------------------------------------------
##
## Script name: 6_cleanNmerge.R
##
## Purpose of script: Bring FT and MO analysis outputs together and clean.
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
## Notes: This script creates figure S3.
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
pckgs <- list("phyloseq", 
              "plyr", "tidyverse", "data.table", # wrangling
              "ape") # change as needed

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

# output from model classification
mo <- read.csv("./Data/Summary/MO_patterns.classified_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()
ft <- read.csv("./Data/Summary/FT_patterns.classified_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()

pks.df <- readRDS(select_newest("./Objects", "MO_pks.df_binned_")) %>% setDT()
mo[pks.df, no.peak := i.no.peak, on = .(ID)]
pks.df <- readRDS(select_newest("./Objects", "FT_pks.df_binned_")) %>% setDT()
ft[pks.df, no.peak := i.no.peak, on = .(ID)]

# # read in database tree
# pb <- readRDS("./Data/Microbial/phyloseq_otu.tax.samp.phylo_2015-18_SINA.rds")
# otus <- t(otu_table(pb))
# # subset OTUs
# sub.otus <- subset(otus, rownames(otus) %in% unique(sapply(str_split(mo$ID, "[.]"), "[[",4)))
# # subset samples
# cleaned.meta <- read.csv('./Data/Microbial/cleaned.mic_meta.df.csv', stringsAsFactors = F)
# samples <- unique(cleaned.meta$seq_name)
# meta <- sample_data(pb)
# sub.samples <- subset(meta, rownames(meta) %in% samples)
# # extract final phyloseq object
# 
# # specify root of tree
# # Function from: john-quensen.com/r/unifrac-and-tree-roots/
# root.outgroup <- function(unrooted.tree){
#   treeDT <- cbind(data.table(unrooted.tree$edge),
#                   data.table(length = unrooted.tree$edge.length))[1:ape::Ntip(unrooted.tree)] %>%
#     cbind(data.table(id = unrooted.tree$tip.label))
#   # Take the longest terminal branch as outgroup
#   new.outgroup <- treeDT[which.max(length)]$id
#   return(new.outgroup)
# }
# 
# outgroup <- root.outgroup(phy_tree(pb))
# rooted.tree <- ape::root(phy_tree(pb), outgroup = outgroup, resolve.root = T)
# 
# physeq <- merge_phyloseq(sub.otus, tax_table(pb), sub.samples, rooted.tree)
# saveRDS(physeq, "./Data/phyloseq_otu.tax.samp.phylo_SINA_chapter2.rds")
# ape::write.tree(phy_tree(physeq), file = "./Data/Microbial/chapter2_rooted_SINA.tree")

# if tree figrue is used in publication cite
# # Letunic I and Bork P (2021) Nucleic Acids Res doi: 10.1093/nar/gkab301 Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation

# all models as list
# ft.ls <- readRDS(select_newest("./Objects","FT_model.ls"))
# mo.ls <- readRDS(select_newest("./Objects","MO_model.ls"))

# 3. Clean data ----------------------------------------------------------------
## Add factor identifiers ------------------------------------------------------
ft[, dataset := "DOM"]

## split ID info
mo[dna.type == "DNA", dataset := "DNA"]
mo[dna.type == "RNA", dataset := "RNA"]

## Add travel time ranges ------------------------------------------------------
# extract mouth travel time
tt.mouth <- read.csv("./Data/Traveltime/travel.time_byyear_atmouth.csv", sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()
tt.mouth[, season := factor(Month, levels = c("6","8","10"), labels = c("Spring","Summer","Autumn"))]
tt.mouth[, c("year","Year") := list(Year, NULL)]
mo[tt.mouth, tt.at.mouth := i.mouth.trav.time_d, on = .(year, season)]
ft[tt.mouth, tt.at.mouth := i.mouth.trav.time_d, on = .(year, season)]

## Add travel time range categories
ft[x.start == -50 & x.end <= tt.at.mouth, tt.range := "soil-mouth"]
ft[x.start >= 0 & x.end <= tt.at.mouth, tt.range := "stream-mouth"]
ft[x.start == -50 & x.end > tt.at.mouth, tt.range := "soil-lake"]
ft[x.start >= 0 & x.end > tt.at.mouth, tt.range := "stream-lake"]

mo[x.start == -50 & x.end <= tt.at.mouth, tt.range := "soil-mouth"]
mo[x.start >= 0 & x.end <= tt.at.mouth, tt.range := "stream-mouth"]
mo[x.start == -50 & x.end > tt.at.mouth, tt.range := "soil-lake"]
mo[x.start >= 0 & x.end > tt.at.mouth, tt.range := "stream-lake"]

# identify one-peaking non-linear models
mo[no.peak > 1 & AU == "non-linear increase", AU := "multimodal increase"]
ft[no.peak > 1 & AU == "non-linear increase",  AU := "multimodal increase"]

mo[no.peak > 1 & AU == "non-linear decrease", AU := "multimodal increase"]
ft[no.peak > 1 & AU == "non-linear decrease", AU := "multimodal increase"]

## rename trends category ------------------------------------------------------
ft[ , sAU := factor(AU, levels = c("n.s.", "increase", "non-linear increase",
                                   "increasing quadratic", "increasing bimodal", "multimodal increase",
                                   "unimodal",
                                   "multimodal decrease",
                                   "decreasing bimodal", "decreasing quadratic",
                                   "non-linear decrease", "decrease"), 
                    labels =c("n.s.", "increase", "non-linear increase",
                              "non-linear increase", "multimodal increase", "multimodal increase",
                              "unimodal",
                              "multimodal decrease",
                              "multimodal decrease", "non-linear decrease",
                              "non-linear decrease", "decrease"))]

mo[ , sAU := factor(AU, levels = c("n.s.", "flat","increase", "non-linear increase",
                                   "increasing quadratic", "increasing bimodal", "multimodal increase",
                                   "unimodal",
                                   "multimodal decrease",
                                   "decreasing bimodal", "decreasing quadratic",
                                   "non-linear decrease", "decrease"), 
                    labels =c("n.s.", "flat","increase", "non-linear increase",
                              "non-linear increase", "multimodal increase", "multimodal increase",
                              "unimodal",
                              "multimodal decrease",
                              "multimodal decrease", "non-linear decrease",
                              "non-linear decrease", "decrease"))]

## Remove a few taxa -----------------------------------------------------------
phy <- readRDS("./Data/Microbial/phyloseq_otu.tax.samp.phylo_SINA_chapter2.rds")
tax <- tax_mat(phy)
tax <- as.data.frame(tax) %>% setDT(keep.rownames = "OTU")
mo <- merge(mo, tax %>% dplyr::select(OTU:domain), by = "OTU", all.x = T)

# number of Archaea and Unclassified
length(unique(mo[domain == "Archaea",]$OTU)) # 55
length(unique(mo[is.na(domain),]$OTU)) #2019
length(unique(mo$OTU)) # out of 6103

# How many of significant models were either one of these categories
nrow(mo[is.na(domain) | domain == "Archaea",][c.ns.s == "sig.",]) # 486 significant models were either archaea or Unclassified
nrow(mo[c.ns.s == "sig.",]) # 1934
(468*100) / 1711 # that's 27%

# Remove
mo <- mo[!(is.na(domain) | domain == "Archaea"),]

## How many of the models actually didn't work? -----------------------------
nrow(mo[sAU == "flat",]) #1782
nrow(mo[sAU == "n.s.",]) #608

mo[(n = .N), by = .(dna.type)]
nrow(mo[(sAU == "flat" | sAU == "n.s.") & dna.type == "DNA",]) #1924
nrow(mo[(sAU == "flat" | sAU == "n.s.") & dna.type == "RNA",]) #466

nrow(ft[sAU == "n.s.",]) #211

nrow(mo); nrow(ft)

# Just remove
mo <- mo[sAU != "n.s.",][sAU != "flat",]
ft <- ft[sAU != "n.s.",]

## Add a minimal spatial pattern column -------------------------------------------------------------------------------
mo[, minAU := factor(sAU,
                     levels = c('increase', "non-linear increase", "multimodal increase",
                                'unimodal',
                                'multimodal decrease', 'non-linear decrease', 'decrease'),
                     labels = c('increase','increase','increase','unimodal',
                                'decrease', 'decrease', 'decrease'))]
ft[, minAU := factor(sAU,
                     levels = c('increase', "non-linear increase", "multimodal increase",
                                'unimodal',
                                'multimodal decrease', 'non-linear decrease', 'decrease'),
                     labels = c('increase','increase','increase','unimodal',
                                'decrease', 'decrease', 'decrease'))]

## Get quantiles of slopes ------------------------------------------------
# histogram, calculate quantiles
ft.quar.pos <- quantile(ft[sAU %in% c('increase', 'non-linear increase') & c.ns.s == 'sig.' & lm.slope > 0, ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))
ft.quar.neg <- quantile(ft[sAU %in% c('decrease', 'non-linear decrease') & c.ns.s == 'sig.' & lm.slope < 0, ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))

(ft.hist <- ggplot(ft[sAU %in% c('increase', 'non-linear increase', 
                                 'non-linear decrease', 'decrease') & c.ns.s == "sig.",],
                   aes(x = lm.slope)) +
    theme_bw() +
    geom_histogram(bins = 100, fill = "white", colour = "gray20") +
    geom_vline(xintercept = ft.quar.pos, colour = "tomato", linetype = 'dashed') +
    geom_vline(xintercept = ft.quar.neg, colour = "royalblue", linetype = 'dashed') +
    labs(x = "z-scaled slope of sign. linear models", y = "Frequency", title = "DOM"))

mo.quar.pos <- quantile(mo[sAU %in% c('increase', 'non-linear increase') & c.ns.s == 'sig.' & 
                             lm.slope > 0 & dna.type == "DNA", ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))
mo.quar.neg <- quantile(mo[sAU %in% c('decrease', 'non-linear decrease') & c.ns.s == 'sig.' & lm.slope < 0 &
                             dna.type == "DNA", ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))


(mo.hist <- ggplot(mo[sAU %in% c('increase', 'non-linear increase', 
                                 'non-linear decrease', 'decrease') & c.ns.s == "sig.",],
                   aes(x = lm.slope)) + 
    theme_bw() +
    geom_histogram(bins = 100, fill = "white", colour = "gray20") +
    geom_vline(xintercept = mo.quar.pos, colour = "tomato", linetype = 'dashed') +
    geom_vline(xintercept = mo.quar.neg, colour = "royalblue", linetype = 'dashed') +
    labs(x = "z-scaled slope of sign. linear models", y = "Frequency", title = "MO"))

# ggsave("./Figures/Analysis/slope_histogram_lm_only_detailed.png",ggarrange(mo.hist, ft.hist, nrow = 2),
#        height = 12, width = 15, units = "cm")

# add bins to main data frame
mo[lm.slope <= mo.quar.neg[2], bin.cat := -3]
mo[lm.slope > mo.quar.neg[2] & lm.slope <= mo.quar.neg[3], bin.cat := -2]
mo[lm.slope > mo.quar.neg[3] & lm.slope <= 0, bin.cat := -1]
mo[lm.slope >= 0 & lm.slope <= mo.quar.pos[2], bin.cat := 1]
mo[lm.slope > mo.quar.pos[2] & lm.slope <= mo.quar.pos[3], bin.cat := 2]
mo[lm.slope > mo.quar.pos[3], bin.cat := 3]
# sanity check
mo[is.na(bin.cat) & !is.na(lm.slope) & dna.type == "DNA",]

ft[lm.slope <= ft.quar.neg[2], bin.cat := -3]
ft[lm.slope > ft.quar.neg[2] & lm.slope <= ft.quar.neg[3], bin.cat := -2]
ft[lm.slope > ft.quar.neg[3] & lm.slope <= 0, bin.cat := -1]
ft[lm.slope >= 0 & lm.slope <= ft.quar.pos[2], bin.cat := 1]
ft[lm.slope > ft.quar.pos[2] & lm.slope <= ft.quar.pos[3], bin.cat := 2]
ft[lm.slope > ft.quar.pos[3] , bin.cat := 3]
# sanity check
ft[is.na(bin.cat) & !is.na(lm.slope),]

mo[, bin.cat := factor(bin.cat, levels = c("-3","-2","-1", "1", "2", "3"),
                       labels = c("Q1","Q2", "Q3","Q4","Q5","Q6"))]
ft[, bin.cat := factor(bin.cat, levels = c("-3","-2","-1", "1", "2", "3"),
                       labels = c("Q1","Q2", "Q3","Q4","Q5","Q6"))]

## Reorder -----------------------------------------------------------------------------------
mo <- mo %>% dplyr::select(dataset, ID, OTU, domain, year:dna.type, tt.at.mouth, tt.range, 
                           model:hump, no.peak, ns.s:c.ns.s, AU, sAU, minAU, bin.cat)
ft <- ft %>% dplyr::select(dataset, ID, MF, year:season, tt.at.mouth, tt.range, 
                           model:hump, no.peak, ns.s:c.ns.s, AU, sAU, minAU, bin.cat)

# 4. Save clean version ---------------------------------------------------------------------------------
write.table(mo, paste0("./Data/Summary/MO_patterns.classified_cleaned_",Sys.Date(), ".csv"), sep = ",",
            row.names = F)
write.table(ft, paste0("./Data/Summary/FT_patterns.classified_cleaned_",Sys.Date(), ".csv"), sep = ",",
            row.names = F)


# Plot examples -----------------------------------------------------------------------

# read in microbial models
model.ls <- readRDS("./Objects/MO_model.ls_binned_2024-02-10.rds")

# Create an example plot for DOM and DNA
# examples were picked from a prior run where random examples of
# DNA and DOM models were individually plotted

# DNA representatives:
# Increase 2015.Summer.DNA.OTU_10979
# Decrease 2015.Spring.DNA.OTU_7568

# create colour scale
cols <- RColorBrewer::brewer.pal(n = 8, 'RdBu')
names(cols) <- 8:1

mo <- c('2015.Spring.DNA.OTU_7568','2015.Summer.DNA.OTU_10979')

m.plots <- lapply(model.ls[names(model.ls) %in% mo], function(x){
  newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
  newd$y.pred <- predict(x$model, newd)
  
  (p <- ggplot()+
      theme_bw() +
      theme(legend.position = "right", panel.grid.minor.y = element_blank(),
            legend.box.margin = margin(2)) +
      scale_x_continuous(minor_breaks = seq(-50,1000, 50)) +
      geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
      geom_point(aes(x = x$binned$x, y = x$binned$y,
                     fill = factor(x$binned$n, levels = 1:8)),
                 size = 5, shape = 21) +
      scale_fill_manual(values = rev(cols[names(cols) %in% x$binned$n]), name = "Sample size") +
      geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
      labs(x = "Flow-weighted water age (d)",
           y = "CSS Reads (z-scaled)"))
  return(p)
})
rm(model.ls)


# Load DOM models
model.ls <- readRDS("./Objects/FT_model.ls_binned_2024-02-11.rds")

# DOM represntatives:
# Increase 2015_Spring_C14H9O9N1S0
# Decrease 2015_Spring_C24H21O13N1S0
ft <- list('2015_Spring_C14H9O9N1S0','2015_Spring_C24H21O13N1S0')

f.plots <- lapply(model.ls[names(model.ls) %in% ft], function(x){
  newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
  newd$y.pred <- predict(x$model, newd)
  
  (p <- ggplot()+
      theme_bw() +
      theme(legend.position = "right", panel.grid.minor.y = element_blank(),
            legend.box.margin = margin(2)) +
      scale_x_continuous(minor_breaks = seq(-50,1000, 50)) +
      geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
      geom_point(aes(x = x$binned$x, y = x$binned$y,
                     fill = factor(x$binned$n, levels = 1:8)),
                 size = 5, shape = 21) +
      scale_fill_manual(values = rev(cols[names(cols) %in% x$binned$n]), name = "Sample size") +
      geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
      labs(x = "Flow-weighted water age (d)",
           y = "Peak intensity (z-scaled)"))
  return(p)
})

# create plot just to extract legend
(leg.p <- ggplot(data.frame(x = 1:8, y = rep(1, times = 8)), aes(x, y)) + 
    theme_bw() +
    theme(legend.box.margin = margin(20)) +
    ylim(0, 1) + 
    geom_point(aes(fill = factor(x)), shape = 21, size = 3) +
    scale_fill_brewer(palette = 'RdBu', direction = -1, name = "Sample size"))

mp <- ggarrange(m.plots[[2]] + guides(fill = 'none'), m.plots[[1]] + guides(fill ='none'), 
                nrow = 1, ncol = 2,legend = 'right', labels = c('a','b'))
mp <- annotate_figure(mp, left = text_grob('DNA', face = 'bold', size = 15, rot = 90))
fp <- ggarrange(f.plots[[1]] + guides(fill = 'none'), f.plots[[2]] + guides(fill ='none'), 
                nrow = 1, ncol = 2,legend = 'right', labels = c('c','d'))
fp <- annotate_figure(fp, left = text_grob('DOM', face = 'bold', size = 15, rot = 90))

(p <- ggarrange(mp, fp, nrow = 2, common.legend = T, 
                legend.grob = get_legend(leg.p), legend = 'right'))

# ggsave('Figures/figS3.pdf', p + theme(plot.background = element_rect(fill = 'white')), 
#        height = 15, width = 20, units = 'cm', dpi = 300)

# Make tables for manuscript ----------------------------------------------------------

# Table S2: Number of samples usde to create spatial models of microbial OTUs and DOM MF per campaign
otu.tab <- select_newest("./Data/Microbial", "201520162017_fin_css_otu99_phantomcor_no.dups_")
otu.tab <- as.matrix(read.csv(
  otu.tab,
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))

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

# make sure we only count those we are using
met.df <- met.df[met.df$row.names %in% row.names(otu.tab),]

met.df[,.(n = .N), by = .(dna_type, year, season)]


# row orders need to match between tax.tab and otu.tab
otu.tab <- otu.tab[order(row.names(otu.tab)),]
knitr::kable(permanova.df, "latex", booktabs = T) %>%
  kable_styling(position = "center", full_width = F) %>%
  add_header_above(c(" " = 2, "PERMANOVA" = 4, "PERMDISP" = 3)) %>%
  collapse_rows(columns = 1, valign = "middle")


# Table S3: Number of overall (n) and unique microbial OTUs and DOM molecular formulae per campaign

# combine tables into one
permanova.df <- rbind(perm.dna, perm.dr)

# replace NAs with ""
permanova.df <- permanova.df %>% replace(is.na(.), "")

knitr::kable(permanova.df, "latex", booktabs = T) %>%
  kable_styling(position = "center", full_width = F) %>%
  add_header_above(c(" " = 2, "PERMANOVA" = 4, "PERMDISP" = 3)) %>%
  collapse_rows(columns = 1, valign = "middle")



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
#   [1] ade4_1.7-23         ape_5.8-1           vegan_2.7-2         permute_0.9-10      kableExtra_1.4.0   
# [6] cowplot_1.2.0       phyloseq_1.54.2     BiocManager_1.30.27 plotly_4.12.0       gridExtra_2.3      
# [11] doMC_1.3.8          mgcv_1.9-4          nlme_3.1-168        rstatix_0.7.3       doParallel_1.0.17  
# [16] iterators_1.0.14    foreach_1.5.2       ggpubr_0.6.3        plyr_1.8.9          lubridate_1.9.5    
# [21] forcats_1.0.1       stringr_1.6.0       dplyr_1.2.0         purrr_1.2.1         readr_2.2.0        
# [26] tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.2       tidyverse_2.0.0     data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] rlang_1.1.7         magrittr_2.0.4      otel_0.2.0          compiler_4.5.2      systemfonts_1.3.1  
# [6] vctrs_0.7.1         reshape2_1.4.5      pkgconfig_2.0.3     crayon_1.5.3        fastmap_1.2.0      
# [11] backports_1.5.0     XVector_0.50.0      labeling_0.4.3      CoprManager_0.5.8   rmarkdown_2.30     
# [16] tzdb_0.5.0          ragg_1.5.0          xfun_0.56           jsonlite_2.0.0      biomformat_1.38.0  
# [21] rhdf5filters_1.22.0 Rhdf5lib_1.32.0     broom_1.0.12        cluster_2.1.8.2     R6_2.6.1           
# [26] stringi_1.8.7       RColorBrewer_1.1-3  car_3.1-5           Rcpp_1.1.1          Seqinfo_1.0.0      
# [31] knitr_1.51          IRanges_2.44.0      Matrix_1.7-4        splines_4.5.2       igraph_2.2.2       
# [36] timechange_0.4.0    tidyselect_1.2.1    rstudioapi_0.18.0   abind_1.4-8         codetools_0.2-20   
# [41] lattice_0.22-9      Biobase_2.70.0      withr_3.0.2         S7_0.2.1            evaluate_1.0.5     
# [46] survival_3.8-6      xml2_1.5.2          Biostrings_2.78.0   pillar_1.11.1       carData_3.0-6      
# [51] stats4_4.5.2        generics_0.1.4      S4Vectors_0.48.0    hms_1.1.4           scales_1.4.0       
# [56] glue_1.8.0          lazyeval_0.2.2      tools_4.5.2         ggsignif_0.6.4      rhdf5_2.54.1       
# [61] grid_4.5.2          Formula_1.2-5       cli_3.6.5           textshaping_1.0.4   viridisLite_0.4.3  
# [66] svglite_2.2.2       gtable_0.3.6        digest_0.6.39       BiocGenerics_0.56.0 htmlwidgets_1.6.4  
# [71] farver_2.1.2        htmltools_0.5.9     multtest_2.66.0     lifecycle_1.0.5     httr_1.4.8         
# [76] MASS_7.3-65  