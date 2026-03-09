## -------------------------------------------------------------------------
##
## Script name: 7_proportion.R
##
## Purpose of script: Calculate proportion of unreactive and reactive
##                    microbial and DOM moieties.
##                    Calculate proportion of spatial trends by year and
##                    season. Prepare individual plots to be combined
##                    into a composite plot (Fig.2) in Inkscape.
##
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
## Notes: Some of the plots are heavy and need a LOT of time to be plotted.
##        This script produced individual parts of the composite plot Fig.2
##        This script also produces Table S3.
##
## -------------------------------------------------------------------------

## Use R project with regular scripts, all paths are relative 

# Server set-up -----------------------------------------------------------
## Working directory is set from where the job is submitted
## Load library path, if on a server
# .libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )

# R-setup -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
pckgs <- list("phyloseq", "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr", "plotly", "ggforce","waffle", "gridExtra", "cowplot", # plotting,
              "kableExtra",
              # "xlsx", # making tables for manuscript (LaTeX) and export as excel (for ISME)
              "doMC", # parallel computing
              "vegan", "ape", "ade4", "rstatix", "mgcv", "picante") # change as needed

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

# remotes::install_github("hrbrmstr/waffle")

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
mo <- read.csv("./Data/Summary/MO_patterns.classified_cleaned_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()
ft <- read.csv("./Data/Summary/FT_patterns.classified_cleaned_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()

# Format data ------------------------------------------------------------
# Define spatial pattern (AU) categories ---------------------------------
mo[, intAU := factor(sAU, levels = c("increase", "non-linear increase", "multimodal increase",
                                     "unimodal",
                                     "multimodal decrease", "non-linear decrease", "decrease"),
                     labels = c("increase", "non-linear increase",'non-linear increase',
                                'unimodal',
                                'non-linear decrease', 'non-linear decrease', 'decrease'))]

ft[, intAU := factor(sAU, levels = c("increase", "non-linear increase", "multimodal increase",
                                     "unimodal",
                                     "multimodal decrease", "non-linear decrease", "decrease"),
                     labels = c("increase", "non-linear increase",'non-linear increase',
                                'unimodal',
                                'non-linear decrease', 'non-linear decrease', 'decrease'))]

# Table S3 --------------------------------------------------------------
## How many MF/OTU in each season / year? -------------------------------
camp.n.ft <- ft[, .(n.camp = .N), by = .(year, season)]
camp.n.mo <- mo[!(dna.type == "RNA"), .(n.camp = .N), by = .(year, season)]

## How many are unique within each season / year -------------------------
# Get all unique MFs
mf.ls <- dlply(ft, .(year, season), function(x){
  unique(x$MF)
})

## Count how many MF/OTU are unique within each campaign -----------------
ft.uni <- data.frame()
names(mf.ls) # order
for(i in 1:length(mf.ls)){
  ft.uni[i,"ID"]<- names(mf.ls)[i]
  ft.uni[i, "uni.n"] <- length(mf.ls[[i]][!(mf.ls[[i]] %in% unique(Reduce(c, mf.ls[-i])))])
}

ft.uni<- ft.uni %>% separate(col = "ID", into = c("year","season")) %>% setDT()

# same for microbes
# Get all unique MFs
mf.ls <- dlply(mo[!(dna.type == "RNA"),], .(year, season), function(x){
  unique(x$OTU)
})

# Count how many are unique within each campaign
mo.uni <- data.frame()
names(mf.ls) # order
for(i in 1:length(mf.ls)){
  mo.uni[i,"ID"]<- names(mf.ls)[i]
  mo.uni[i, "uni.n"] <- length(mf.ls[[i]][!(mf.ls[[i]] %in% unique(Reduce(c, mf.ls[-i])))])
}

mo.uni<- mo.uni %>% separate(col = "ID", into = c("year","season")) %>% setDT()

# format data
ft.uni[, year := as.numeric(year)]
mo.uni[, year := as.numeric(year)]

# merge together
camp.n.ft <- camp.n.ft[ft.uni, , on = .(year, season)]
camp.n.mo <- camp.n.mo[mo.uni, , on = .(year, season)]

camp.n <- bind_rows(camp.n.ft[, dataset := "DOM"],
                    camp.n.mo[, dataset := "MO"]) %>% setDT()

# format for table
temp <- dcast(camp.n, year + season ~ dataset, value.var = c("n.camp","uni.n"))
temp <- temp %>% dplyr::select(year:season, n.camp_MO, uni.n_MO, n.camp_DOM, uni.n_DOM)

# change column names
colnames(temp) <- c("Year","Season","n","unique","n","unique")

# format/clean output in LaTeX
knitr::kable(temp, "latex", booktabs = T) %>%
  kable_styling(position = "center", full_width = F) %>%
  add_header_above(c(" " = 2, "Microbial" = 2, "DOM" = 2))

rm(camp.n.ft, camp.n.mo, camp.n, temp, ft.uni, mo.uni, mf.ls)

# Proportion of classified trends -----------------------------------------------------------------
## overall percentage unreactive vs. reactive per year and season ----------------------------------------------------
ft[, n := .N, by = .(year, season)]
nsft <- ft[,.(n = unique(n),
              n.ns.s = .N), by = .(c.ns.s, year, season)]
nsft[, perc := round((n.ns.s * 100) / n, 0)]

mo[dna.type != "DNA", n := .N, by = .(year, season)]
nsmo <- mo[dna.type != "DNA",.(n = unique(n),
                               n.ns.s = .N), by = .(c.ns.s, year, season)]
nsmo[, perc := round((n.ns.s * 100) / n, 0)]
nsmo <- nsmo[order(c.ns.s),]

mo[dna.type != "RNA", n := .N, by = .(year, season)]
nsmo <- mo[dna.type != "RNA",.(n = unique(n),
                               n.ns.s = .N), by = .(c.ns.s, year, season)]
nsmo[, perc := round((n.ns.s * 100) / n, 0)]
nsmo <- nsmo[order(c.ns.s),]


nsmo <- mo[dna.type != "RNA",.(n = unique(n),
                               n.ns.s = .N), by = .(c.ns.s, year, season)]
nsmo[, perc := round((n.ns.s * 100) / n, 0)]
nsmo <- nsmo[order(c.ns.s),]

nsdf.year <- bind_rows(nsft[, dataset := "DOM"], nsmo[, dataset := "MO"])
# san check
nsmo[, .(sum = sum(perc)), by = .(year, season)]

waf.ft <- ggplot(nsdf.year[dataset == "DOM",]) +
  facet_grid(season~year) +
  geom_waffle(aes(fill = c.ns.s, values = perc), n_rows = 4, size = 0.33, colour = 'white', na.rm = T) +
  scale_fill_manual(name = NULL, values = c("gray70", "gray10"),
                    labels = c("Unreactive", "Reactive")) +
  coord_equal() +
  theme_void() +
  labs(title = "", xlab = "") + #title = "DOM", xlab = "1 square = 1 %"
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 15), #20
        legend.justification = c(0,0.85),
        axis.title.x = element_text(hjust = 0.9)) #plot.margin = unit(c(0, 0, 0, 0), "cm")

waf.mo <- ggplot(nsdf.year[dataset == "MO",], aes(fill = c.ns.s, values = perc)) +
  facet_grid(season~year) +
  waffle::geom_waffle(n_rows = 4, size = 0.33, colour = 'white', na.rm = T) +
  scale_fill_manual(name = NULL, values = c("gray70", "gray10"),
                    labels = c("Unreactive", "Reactive")) +
  labs(title = "", xlab = "") + #title = "MO",
  coord_equal() +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 15), #20
        legend.justification = c(0,0.85),
        axis.title.x = element_text(hjust = 0.9)) #plot.margin = unit(c(0, 0, 0, 0), "cm")

# ggsave("./Figures/Analysis/waffle_perc_byyearseason_mic.png", waf.mo,
#        height = 10, width = 15, units = 'cm')
# ggsave("./Figures/Analysis/waffle_perc_byyearseason_ft.png", waf.ft,
#        height = 10, width = 15, units = 'cm')

## overall percentage unreactive vs. reactive ---------------------------
# get total N
ft[, n := .N]
# get number of unreactive models
nsft <- ft[,.(n = unique(n),
              n.ns.s = .N), by = .(c.ns.s)]
# calculate percentage
nsft[, perc := round((n.ns.s * 100) / n, 0)]

# get number of DNA models
mo[dna.type != "DNA", n := .N]
# get number of unreactive models
nsmo <- mo[dna.type != "DNA",.(n = unique(n),
                               n.ns.s = .N), by = .(c.ns.s)]
# calculate percentage
nsmo[, perc := round((n.ns.s * 100) / n, 0)]
nsmo <- nsmo[order(c.ns.s),]

# get number of RNA models
mo[dna.type != "RNA", n := .N]
# get number of unreactive models
nsmo <- mo[dna.type != "RNA",.(n = unique(n),
                               n.ns.s = .N), by = .(c.ns.s)]
# calculate percentage
nsmo[, perc := round((n.ns.s * 100) / n, 0)]
nsmo <- nsmo[order(c.ns.s),]

# merge DOM and MO together
nsdf <- bind_rows(nsft[, dataset := "DOM"], nsmo[, dataset := "MO"])

waf.ft <- ggplot(nsdf[dataset == "DOM",]) +
  #facet_grid(season~year) +
  geom_waffle(aes(fill = c.ns.s, values = perc), n_rows = 4, size = 0.33, colour = 'white', na.rm = T) +
  scale_fill_manual(name = NULL, values = c("gray70", "gray10"),
                    labels = c("Unreactive", "Reactive")) +
  coord_equal() +
  theme_void() +
  labs(title = "", xlab = "") + #title = "DOM", xlab = "1 square = 1 %"
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 15), #20
        legend.justification = c(0,0.85),
        axis.title.x = element_text(hjust = 0.9)) #plot.margin = unit(c(0, 0, 0, 0), "cm")

waf.mo <- ggplot(nsdf[dataset == "MO",], aes(fill = c.ns.s, values = perc)) +
  #facet_grid(season~year) +
  waffle::geom_waffle(n_rows = 4, size = 0.33, colour = 'white', na.rm = T) +
  scale_fill_manual(name = NULL, values = c("gray70", "gray10"),
                    labels = c("Unreactive", "Reactive")) +
  labs(title = "", xlab = "") + #title = "MO",
  coord_equal() +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 15), #20
        legend.justification = c(0,0.85),
        axis.title.x = element_text(hjust = 0.9)) #plot.margin = unit(c(0, 0, 0, 0), "cm")

wafs <- plot_grid(waf.mo + theme(legend.position = "none"),
                  waf.ft + theme(legend.position = "none"),
                  align = "vh", nrow = 1, hjust = -1)

#ggsave("./Figures/fig2_moft.waf.png", wafs, height = 5, width = 10, units = 'cm')

waf.fin <- plot_grid(wafs, get_legend(waf.ft), nrow = 1, ncol = 2, align = "v", axis = "tb", rel_widths = c(1,0.1))

# these plots will be part of composite plot. Save as .eps and merge the plots in Inkscape.

# save_plot("./Figures/fig2_mo.waf.eps", 
#           waf.mo + theme(legend.position = "none",
#                          panel.background = element_rect(fill = "transparent",
#                                                          colour = NA_character_), # necessary to avoid drawing panel outline
#                          panel.grid.major = element_blank(), # get rid of major grid
#                          panel.grid.minor = element_blank(), # get rid of minor grid
#                          plot.background = element_rect(fill = "transparent",
#                                                         colour = NA_character_), # necessary to avoid drawing plot outline
#                          legend.background = element_rect(fill = "transparent"),
#                          legend.box.background = element_rect(fill = "transparent"),
#                          legend.key = element_rect(fill = "transparent")), dpi = 300,
#           base_height = 2, base_width = 5)
# 
# save_plot("./Figures/fig2_ft.waf.eps", 
#           waf.ft + theme(legend.position = "none",
#                          panel.background = element_rect(fill = "transparent",
#                                                          colour = NA_character_), # necessary to avoid drawing panel outline
#                          panel.grid.major = element_blank(), # get rid of major grid
#                          panel.grid.minor = element_blank(), # get rid of minor grid
#                          plot.background = element_rect(fill = "transparent",
#                                                         colour = NA_character_), # necessary to avoid drawing plot outline
#                          legend.background = element_rect(fill = "transparent"),
#                          legend.box.background = element_rect(fill = "transparent"),
#                          legend.key = element_rect(fill = "transparent")), dpi = 300,
#           base_height = 2, base_width = 5)
# 
# save_plot("./Figures/fig2_waf.legend.eps", as_ggplot(get_legend(waf.ft + theme(legend.position = "top",
#                                                                                             plot.background = element_rect(fill = "transparent",
#                                                                                                                            colour = NA_character_), # necessary to avoid drawing plot outline
#                                                                                             #legend.background = element_rect(fill = "transparent"),
#                                                                                             #legend.box.background = element_rect(fill = "transparent"),
#                                                                                             #legend.key = element_rect(fill = "transparent")), 
# ))), dpi = 300, base_height = 2, base_width = 2)
# 
# ggsave("./Figures/fig2_mo.waf.png", waf.mo + guides(fill = "none"), height = 5, width = 10, units = "cm", dpi = 300)
# ggsave("./Figures/fig2_ft.waf.png", waf.ft + guides(fill = "none"), height = 5, width = 10, units = "cm", dpi = 300)



# Proportion of spatial patterns --------------------------------------------------
# get number of models by year and season (total n)
ft[c.ns.s != "n.s.", n.camp := .N, by = .(year, season)]
mo[c.ns.s != "n.s." & sAU != "n.s." & dna.type != "RNA", n.camp := .N, by = .(year, season)]

# get number models in each spatial pattern
cft <- ft[c.ns.s != "n.s.",.(n = unique(n.camp),
                             n.camp = .N), by = .(intAU, year, season)]
# calculate percentage by spatial pattern
cft[, perc := (n.camp * 100) / n]

# sanity check
cft[, .(sum = sum(perc, na.rm = T)), by = .(year, season)] # should be all 100%

# Do the same for microbes
# get number of models by year and season (total n)
cmo <- mo[c.ns.s != "n.s." & sAU != "n.s." & dna.type != "RNA",.(n = unique(n.camp),
                                                                 n.camp = .N), by = .(intAU, year, season)]
# get proportion
cmo[, perc := (n.camp * 100) / n]

# sanity check
cmo[, .(sum = sum(perc, na.rm = T)), by = .(year, season)] # should be all 100%

# Make an overall pattern proportion plot
perc.camp <- rbind(cft[, dataset := "DOM"],
                   cmo[,dataset := "MO"])

perc.camp <- perc.camp[!is.na(perc.camp$n),]
# set season as factor
perc.camp[, season := factor(season, levels = c("Summer","Spring"))]

# set colour palette
library(RColorBrewer)
col.vec <- brewer.pal(5, "RdBu")
# give colour vector names of pattern types
names(col.vec) <- levels(ft$intAU)

(bar.ft <- ggplot(perc.camp[dataset == "DOM",], aes(y = season, x = perc)) +
    theme_pubr() +
    facet_wrap(~year, scales = "free_y", ncol = 1,  strip.position = "left") +
    scale_fill_manual(values = rev(col.vec[names(col.vec) %in% 
                                             unique(perc.camp[dataset == "DOM",]$intAU)]), name = "Spatial patterns",
                      guide = guide_legend(reverse = TRUE)) +
    geom_col(aes(fill = forcats::fct_rev(intAU)), colour = "gray20") +
    scale_x_continuous(position = "top") +
    labs(x = "Percentage within significant models (%)", y = NULL) +
    #coord_flip() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 10), #angle = 90, vjust = 0.5
          axis.title = element_text(size = 12),
          legend.justification = c(0,0.5),
          strip.text = element_text(colour = "gray20", size = 14),
          legend.position = "right"))

(bar.mo <- ggplot(perc.camp[dataset == "MO",], aes(y = season, x = perc)) +
    theme_pubr() +
    facet_wrap(~year, scales = "free_y",ncol = 1,  strip.position = "left") +
    scale_fill_manual(values = rev(col.vec[names(col.vec) %in% unique(perc.camp[dataset == "MO",]$intAU)]), name = "Spatial patterns",
                      guide = guide_legend(reverse = TRUE)) +
    geom_col(aes(fill = forcats::fct_rev(intAU)), colour = "gray20") +
    scale_x_continuous(position = "top") +
    labs(x = "Percentage within significant models (%)", y = NULL) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 10), #, angle = 90, vjust = 0.5
          axis.title = element_text(size = 12),
          strip.text = element_text(colour = "gray20", size = 14),
          legend.position = "bottom"))

bar.for.leg <-ggplot(perc.camp, aes(y = season, x = perc)) +
  theme_pubr() +
  facet_wrap(~year, scales = "free_y",ncol = 1,  strip.position = "left") +
  scale_fill_manual(values = rev(col.vec[names(col.vec) %in% unique(perc.camp$intAU)]), name = "Spatial patterns",
                    guide = guide_legend(reverse = TRUE)) +
  geom_col(aes(fill = forcats::fct_rev(intAU)), colour = "gray20") +
  scale_x_continuous(position = "top") +
  labs(x = "Percentage within significant models (%)", y = NULL) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10), #, angle = 90, vjust = 0.5
        axis.title = element_text(size = 12),
        strip.text = element_text(colour = "gray20", size = 14),
        legend.position = "bottom")

# Again, these will be part of composite plot. Save as eps and make collage in Inkscape.

# save_plot("./Figures/fig2_ft.bars.eps", 
#           bar.ft + theme(legend.position = "none",
#                          panel.background = element_rect(fill = "transparent",
#                                                          colour = NA_character_), # necessary to avoid drawing panel outline
#                          panel.grid.major = element_blank(), # get rid of major grid
#                          panel.grid.minor = element_blank(), # get rid of minor grid
#                          plot.background = element_rect(fill = "transparent",
#                                                         colour = NA_character_), # necessary to avoid drawing plot outline
#                          legend.background = element_rect(fill = "transparent"),
#                          legend.box.background = element_rect(fill = "transparent"),
#                          legend.key = element_rect(fill = "transparent")), dpi = 300,
#           base_height = 2.1, base_width = 5)
# 
# save_plot("./Figures/fig2_mo.bars.eps", 
#           bar.mo + theme(legend.position = "none",
#                          panel.background = element_rect(fill = "transparent",
#                                                          colour = NA_character_), # necessary to avoid drawing panel outline
#                          panel.grid.major = element_blank(), # get rid of major grid
#                          panel.grid.minor = element_blank(), # get rid of minor grid
#                          plot.background = element_rect(fill = "transparent",
#                                                         colour = NA_character_), # necessary to avoid drawing plot outline
#                          legend.background = element_rect(fill = "transparent"),
#                          legend.box.background = element_rect(fill = "transparent"),
#                          legend.key = element_rect(fill = "transparent")), dpi = 300,
#           base_height = 2.1, base_width = 5)
# 
# ggsave("./Figures/fig2_mo.bars.png", bar.mo + guides(fill = "none"),
#        height = 5, width = 12, units = 'cm', dpi = 300)
# ggsave("./Figures/fig2_mo.bars.png", bar.ft + guides(fill = "none"),
#        height = 5, width = 12, units = 'cm', dpi = 300)
# ggsave("./Figures/fig2_bars.legend.png", as_ggplot(get_legend(bar.ft + theme(legend.position = "top",
#                                                                                                         legend.text = element_text(size = 8),
#                                                                                                         legend.title = element_text(size = 10),
#                                                                                                         legend.spacing.x = unit(0.45, 'cm')))),
#        height = 5, width = 19, units = 'cm', dpi = 300)
# 
# # save legend
# save_plot(
#   "./Figures/fig2_bars.legend.eps",
#   as_ggplot(get_legend(
#     bar.for.leg + #guides(fill=guide_legend(override.aes=list(colour =NA), reverse = T))+
#       theme(
#         plot.background = element_rect(fill = "transparent",
#                                        colour = NA_character_),
#         #legend.key = element_rect(colour = "transparent"),
#         # necessary to avoid drawing plot outline
#         legend.background = element_rect(fill = "transparent"),
#         #legend.box.background = element_rect(fill = "transparent"),
#         #legend.key = element_rect(fill = "transparent"))
#       ))
#   ),
#   dpi = 300,
#   base_height = 3,
#   base_width = 8
# )

# A general idea how they will together

p <- plot_grid(waf.mo + theme(legend.position = "none"),
               waf.ft + theme(legend.position = "none"),
               bar.mo + theme(legend.position = "none"),
               bar.ft + theme(legend.position = "none"),
               align = "v", axis = "rl", nrow = 2, ncol = 2, rel_heights = c(0.2,0.8),
               labels = c("a","","b",""))

leg <- plot_grid(get_legend(waf.ft), get_legend(bar.for.leg), ncol = 1)

waf.bar <- plot_grid(p, leg,rel_widths = c(0.75, 0.25), align = "hv", axis = "tbrl")

# save_plot("./Figures/fig2_mo.ft.bars.png", p, dpi = 300,
#           base_height = 4, base_width = 9)
# save_plot("./Figures/fig2_mo.ft.bars.legend.png", leg, dpi = 300,
#           base_height = 4, base_width = 11)


# Pattern plot ---------------------------------------------------------------

# Let's calculate how many models are there for each range
ft[, .(n = .N), by = .(tt.range)]
mo[, .(n = .N), by = .(tt.range)]

# load molecular formula models
ft.model.ls <- readRDS(select_newest("./Objects", "FT_model.ls_binned_"))

ft[c.ns.s == "n.s.", intAU := "n.s."]
ft[, intAU := factor(intAU, levels = c("n.s.",
                                       "increase",
                                       "non-linear increase",
                                       "unimodal",
                                       "non-linear decrease",
                                       "decrease"))]

# loop through each pattern group
groups <- as.vector(levels(ft$intAU)) #[-1]
trav.range <- as.vector(unique(ft$tt.range))

# get the prediction lines for each model
# get 100 examples for each type
set.seed(3)
ft.pred <- list()
for(i in 1:length(groups)){
  for(j in 1:length(trav.range)){
    # get all IDs from the same group
    if(i == 1){
      temp <- ft[c.ns.s == "n.s.",]$ID
      temp <- sample(temp, 100)
    } else {
      temp <- ft[c.ns.s != "n.s." & intAU == groups[i] & tt.range == trav.range[j],]$ID
      if(length(temp) > 100L){
        temp <- sample(temp, 100)
      }
    }
    
    
    # if(length(temp) > 160){
    #   temp <- sample(temp, 160)
    # }
    
    pred.df <- ldply(ft.model.ls[names(ft.model.ls) %in% temp], function(x){
      newd <- data.frame(x =seq(min(x$raw$x, na.rm = T),max(x$raw$x, na.rm = T), by = 1))
      newd$y.pred <- predict(x$model, newd)
      newd$group <- groups[i]
      newd$range <- trav.range[j]
      return(newd)
    })
    
    ft.pred[[paste(groups[i], trav.range[j])]] <- pred.df
  }}

# merge all together
ft.pred <- bind_rows(ft.pred) %>% setDT()

ft.pred[ , group := factor(group, levels = c("n.s.","increase",
                                             "non-linear increase",
                                             "unimodal",
                                             "non-linear decrease",
                                             "decrease"))]

# Same for microbes
# Load microbial models
mo.model.ls <- readRDS(select_newest("./Objects", "MO_model.ls_binned_"))

mo[c.ns.s == "n.s.", intAU := "n.s."]
mo[, intAU := factor(intAU, levels = c("n.s.",
                                             "increase",
                                             "non-linear increase",
                                             "unimodal",
                                             "non-linear decrease",
                                             "decrease"))]
# go through each pattern group
groups <- as.vector(levels(mo$intAU)) #[-1]

# get the prediction lines for each model
# get 100 examples for each type
set.seed(3)
mo.pred <- list()
for(i in 1:length(groups)){
  for(j in 1:length(trav.range)){
    # get all IDs from the same group
    if(i == 1){
      temp <- mo[c.ns.s == "n.s." & sAU != "n.s." & dna.type != "RNA" & tt.range == trav.range[j],]$ID
      if(length(temp) > 100L){
        temp <- sample(temp, 100)
      }
    } else {
      temp <- mo[c.ns.s != "n.s." & intAU == groups[i] & dna.type != "RNA" & tt.range == trav.range[j],]$ID
      if(length(temp) > 100L){
        temp <- sample(temp, 100)
      }
    }
    # if(length(temp) > 160){
    #   temp <- sample(temp, 160)
    # }
    
    pred.df <- ldply(mo.model.ls[names(mo.model.ls) %in% temp], function(x){
      newd <- data.frame(x =seq(min(x$raw$x, na.rm = T),max(x$raw$x, na.rm = T), by = 1))
      pred <- try(predict(x$model, newd), silent = T)
      if (!inherits(pred, "try-error")) { 
        newd$y.pred <- pred
      } else {
        newd$y.pred <- as.numeric(NA)
      }
      newd$group <- groups[i]
      newd$range <- trav.range[j]
      return(newd)
    })
    
    mo.pred[[paste(groups[i], trav.range[j])]] <- pred.df
  }}

# merge all togeher
mo.pred <- bind_rows(mo.pred) %>% setDT()

# set group as factor
mo.pred[ , group := factor(group, levels = c("n.s.",
                                             "increase",
                                             "non-linear increase",
                                             "unimodal",
                                             "non-linear decrease",
                                             "decrease"))]
# get DNA models only
mo.pred <- mo.pred[grep("DNA", .id),]


# Different smoothing methods depending on facet
# Set smoothing methods for each plot panel for the different pattern types
meths <- c("loess","lm","loess","loess",
           #"loess","loess",
           "loess","lm")

# Smoothing function with different behaviour in the different plot panels
custom.smooth <- function(formula,data,...){
  meth <- eval(parse(text=meths[unique(data$PANEL)]))
  x <- match.call()
  x[[1]] <- meth
  eval.parent(x)
}

# set labeller for clean plotting
labeller <- c(`increase` = "Increase",
              `non-linear increase` = "Non-linear\nincrease",
              #`multimodal increase` = "Multimodal\nincrease",
              `unimodal` = "Unimodal",
              #`multimodal decrease` = "Multimodal\ndecrease",
              `non-linear decrease` = "Non-linear\ndecrease",
              `decrease` = "Decrease",
              `n.s.` = "Non-sign.")

# set colour vector
col.vec <- c("#555555","#B2182B","#D6604D",
             "gray70",
             "#4393C3","#2166AC")

# Plot MO patterns

# ! warning ! Very heavy, might take a while to plot.
# Save plot to object but not plot it in R Studio.
# We will look at the plot once it is saved on the PC.

mo.pat <- ggplot(mo.pred)+
  theme_bw() +
  geom_line(aes(x, y.pred, colour = group, group = .id), alpha = 0.5) +
  scale_y_continuous(breaks = seq(-1, 4, 1), limit = c(0,4), position = "left") +
  scale_x_continuous(breaks = c(0,1000)) +
  scale_colour_manual(values = col.vec, name = "Spatial patterns") + #"#FFEE99"
  facet_grid(.~ group, scales = "free_y", labeller = as_labeller(labeller)) + #switch = "y" , range
  labs(x = "Travel time (d)", y = "Reads\n(z-scaled)") +
theme(strip.background = element_rect(fill = "white", colour = "white"), #element_rect(colour = "gray20")
      strip.text.x = element_text(size = 10, colour = "black"), #poster 14
      axis.title.x = element_text(size = 12, colour = "white"),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10, colour = "white"),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.ticks.x = element_line(colour = "white"),
      panel.grid = element_blank())
# plot smooth lines
mo.pat <- mo.pat + geom_smooth(aes(x = x, y = y.pred, linetype = range), method="custom.smooth", 
                               se = F, colour = "gray20", linewidth = 0.5) +
  scale_linetype_manual(values = c('solid','dotted','dashed','twodash'), name = "Ecosystem \nrange") +
  guides(colour= 'none')

# # save for composite plot in Inkscape.
# save_plot("./Figures/fig2_ecosys.legend.eps", as_ggplot(get_legend(mo.pat + theme(legend.position = "right",
#                                                                                             plot.background = element_rect(fill = "transparent",
#                                                                                                                            colour = NA_character_), # necessary to avoid drawing plot outline
#                                                                                             #legend.background = element_rect(fill = "transparent"),
#                                                                                             #legend.box.background = element_rect(fill = "transparent"),
#                                                                                             #legend.key = element_rect(fill = "transparent")), 
# ))), dpi = 300, base_height = 3, base_width = 3)


# Plot MF patterns

# ! warning ! Very heavy, might take a while to plot.
# Save plot to object but not plot it in R Studio.
# We will look at the plot once it is saved on the PC.

ft.pat <- ggplot(ft.pred)+
  theme_bw() +
  #theme(legend.position = "none") +
  #geom_hline(yintercept = 0, colour = "grey20") +
  geom_line(aes(x, y.pred, colour = group, group = .id), alpha = 0.5) +
  scale_y_continuous(breaks = seq(-1, 4, 1), limit = c(0,4), position = "left") +
  scale_x_continuous(breaks = c(0,1000)) +
  scale_colour_manual(values = col.vec, name = "Spatial patterns") + #"#FFEE99"
  facet_grid(. ~ group, scales = "free_y", labeller = as_labeller(labeller)) + #switch = "y"
  labs(x = "Travel time (d)", y = "Peak intensity\n(z-scaled)") +
theme(strip.background = element_rect(fill = "white", colour = "white"), #element_rect(colour = "gray20")
      strip.text.x = element_text(size = 10, colour = "black"), #poster 14
      axis.title.x = element_text(size = 12, colour = "white"),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10, colour = "white"),
      axis.text.y = element_text(size = 10, colour = "black"),
      axis.ticks.x = element_line(colour = "white"),
      panel.grid = element_blank())

# plot smooth lines
ft.pat <- ft.pat + geom_smooth(aes(x = x, y = y.pred, linetype = range), method="custom.smooth", se = F, colour = "gray20",
                               linewidth = 0.5) +
  scale_linetype_manual(values = c('solid','dotted','dashed','twodash'), name = "Ecosystem \nrange") +
  guides(colour= 'none', linetype = "none")


# ggsave("./Figures/fig2_mo.patterns.eps",
#        mo.pat + theme(
#          panel.grid.major = element_blank(), # get rid of major grid
#          panel.grid.minor = element_blank(), # get rid of minor grid
#          plot.background = element_rect(fill = "transparent",
#                                         colour = NA_character_), # necessary to avoid drawing plot outline
#          legend.background = element_rect(fill = "transparent", colour = NA_character_),
#          legend.key = element_rect(fill = "transparent")) + guides(linetype = "none"),
#        dpi = 500, height = 8, width = 23, unit = 'cm',
#        device=cairo_ps, fallback_resolution = 600)
# 
# ggsave("./Figures/fig2_ft.patterns.eps", 
#        ft.pat +
#          theme(
#            panel.grid.major = element_blank(), # get rid of major grid
#            panel.grid.minor = element_blank(), # get rid of minor grid
#            plot.background = element_rect(fill = "transparent",
#                                           colour = NA_character_), # necessary to avoid drawing plot outline
#            legend.background = element_rect(fill = "transparent", colour =  NA_character_),
#            legend.key = element_rect(fill = "transparent")), dpi = 500, height = 8, width = 23, unit = 'cm',
#        device=cairo_ps, fallback_resolution = 600)



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
#   [1] RColorBrewer_1.1-3  picante_1.8.2       mgcv_1.9-4          nlme_3.1-168        rstatix_0.7.3      
# [6] ade4_1.7-23         ape_5.8-1           vegan_2.7-2         permute_0.9-10      doMC_1.3.8         
# [11] iterators_1.0.14    foreach_1.5.2       kableExtra_1.4.0    cowplot_1.2.0       gridExtra_2.3      
# [16] waffle_1.0.2        ggforce_0.5.0       plotly_4.12.0       ggpubr_0.6.3        data.table_1.18.2.1
# [21] lubridate_1.9.5     forcats_1.0.1       stringr_1.6.0       dplyr_1.2.0         purrr_1.2.1        
# [26] readr_2.2.0         tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.2       tidyverse_2.0.0    
# [31] plyr_1.8.9          phyloseq_1.54.2    
# 
# loaded via a namespace (and not attached):
#   [1] rlang_1.1.7         magrittr_2.0.4      otel_0.2.0          compiler_4.5.2      systemfonts_1.3.1  
# [6] vctrs_0.7.1         reshape2_1.4.5      pkgconfig_2.0.3     crayon_1.5.3        fastmap_1.2.0      
# [11] backports_1.5.0     XVector_0.50.0      labeling_0.4.3      rmarkdown_2.30      CoprManager_0.5.8  
# [16] tzdb_0.5.0          xfun_0.56           jsonlite_2.0.0      biomformat_1.38.0   rhdf5filters_1.22.0
# [21] Rhdf5lib_1.32.0     tweenr_2.0.3        broom_1.0.12        cluster_2.1.8.2     R6_2.6.1           
# [26] stringi_1.8.7       car_3.1-5           extrafontdb_1.1     Rcpp_1.1.1          Seqinfo_1.0.0      
# [31] knitr_1.51          extrafont_0.20      IRanges_2.44.0      Matrix_1.7-4        splines_4.5.2      
# [36] igraph_2.2.2        timechange_0.4.0    tidyselect_1.2.1    rstudioapi_0.18.0   abind_1.4-8        
# [41] codetools_0.2-20    curl_7.0.0          lattice_0.22-9      Biobase_2.70.0      withr_3.0.2        
# [46] S7_0.2.1            evaluate_1.0.5      survival_3.8-6      polyclip_1.10-7     xml2_1.5.2         
# [51] Biostrings_2.78.0   pillar_1.11.1       carData_3.0-6       DT_0.34.0           stats4_4.5.2       
# [56] generics_0.1.4      S4Vectors_0.48.0    hms_1.1.4           scales_1.4.0        glue_1.8.0         
# [61] lazyeval_0.2.2      tools_4.5.2         ggsignif_0.6.4      rhdf5_2.54.1        grid_4.5.2         
# [66] Rttf2pt1_1.3.14     Formula_1.2-5       cli_3.6.5           textshaping_1.0.4   viridisLite_0.4.3  
# [71] svglite_2.2.2       gtable_0.3.6        digest_0.6.39       BiocGenerics_0.56.0 htmlwidgets_1.6.4  
# [76] farver_2.1.2        htmltools_0.5.9     multtest_2.66.0     lifecycle_1.0.5     httr_1.4.8         
# [81] MASS_7.3-65  