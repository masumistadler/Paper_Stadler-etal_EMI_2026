## -------------------------------------------------------------------------
##
## Script name: 3_travel.time_2015-2016.R
##
## Purpose of script: We use the path identified in the prior script to actually
##                    calculate flow-weighted travel time or in this case
##                    flow-weighted water age since we included WRT of lentic systems.
##                    
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
## Notes:  This script was run on a high-performance computer 
##         and it is not recommended to try on a personal machine.
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
              "doParallel","doMC","foreach", "doSNOW"
) # change as needed

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

## Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Read in data -----------------------------------------------------------

for(y in 2015:2016) {
  for (mon in c("6", "8")){
    # Read in data
    streamnet <-
      readRDS(paste0("./Objects/flacc3000_streamnet_wWRT_", y, "-", mon, "_final.rds"))
    
    m <- streamnet
    m <- m[!duplicated(m$pointid), ]
    m <- m[water.body == "Fluvial", flag_main.chan := 1]
    
    ## Parallel environment ---------------------------------------------------
    ## Server version
    cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    registerDoMC(cores = cores)

    rm(streamnet)
    
    # Recode Yves' loop in JMP to R -----------------------------------------------------------------------------------------
    
    # Our data is structured as follows
    # next.pixel = indicates which pointid is the next downstream pixel
    # flag_source = indicates whether current pixel is source
    # flag_conf = indicates whether pixel is a confluence
    # flag_end = indicates whether pixel is just before a confluence
    # flag_mouth = indicates whether pixel is the mouth of the river
    # flag_lake.in = indicates whether pixel is the main inlet of a lake/reservoir
    # flag_lake.trib = indicates whether pixel is entering a lake/reservoir
    # flag_lake.out = indicates whether pixel is a outlet of a lake/reservoir
    # velocity_ms.wb = velocity in metres/sec, considering water bodies (reservoir, riverine lakes)
    # wrt_min.wb = water residence time in pixel, considering water bodies (reservoir, riverine lakes)
    
    # First loop ---------------------------------------------------------------------------------
    # Track first stream orders down until a confluence
    # create empty data frame to fill in
    cumtime <- data.frame()
    
    # extract all pointids that are sources
    sources <- m[m$flag_source == 1 & flag_main.chan == 1, ]$pointid
    
    # sum along the path of first order streams
    temp <- foreach(i = 1:length(sources), .combine = rbind) %dopar% {
      
      if (m[m$pointid == sources[i], "flag_source"] == 1) {
        s <- 1
        k <- sources[i]
        cumtime <- data.frame()
        cumtime[s, "pointid"] <- k
        cumtime[s, "travel.time"] <- m[m$pointid == k, "wrt_min.wb"]
        nx.px <- m[m$pointid == k, ]$next.pixel

        # calculate travel time within first order stream
        while (m[m$pointid == nx.px, "flag_conf"] != 1) {
          # while next pixel is not a confluence, do...
          s <- s + 1
          # calculate cumulative travel time by...
          cumtime[s, "pointid"] <- m[m$pointid == k, ]$next.pixel
          cumtime[s, "travel.time"] <-
            cumtime[cumtime$pointid == k, ]$travel.time + # Time A
            (m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$wrt_min.wb * # Time B
               # adding wrt of next pixel to the weighted wrt of the current pixel
               (1 - (m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
                  m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$dis_m3s))
          
          # Calculation:
          # Time A + Time B * (1 - (disB - disA)/disB)
          
          # once this is done, overwrite k with the next pixel
          k <-
            m[m$pointid == k, ]$next.pixel # so the loop jumps to the next pixel
          nx.px <- m[m$pointid == k,]$next.pixel
          }
        }
        cumtime <- cumtime[!is.na(cumtime$travel.time), ]

        return(cumtime)
      }
    
    temp <- temp %>% distinct()
    cumtime <- temp
    saveRDS(cumtime,
            paste0("./Objects/temp_strahord1_", y, "-", mon, ".rds"))
    print("Finished stream orders 1")

    # merge with main streamnetwork file
    temp <- m %>% dplyr::select(pointid)
    cumtime <- merge(temp, cumtime, all = TRUE)
    
    # 2nd loop ------------------------------------------------------------------
    # takes care of confluences and calculates the time in reaches after confluences
    
    # get the maximum stream order in stream network
    max.order <- max(m$strah_ord)
    
    for (s in 2:max.order) {
      print(paste0("Starting loop at stream order ", s))
      #:max.order
      # get all the confluences of the given order
      # order by flow length
      m <- m[order(fllength_dec),]
      conf <- m[m$strah_ord == s & m$flag_conf == 1 & m$flag_main.chan == 1, ]$pointid
      
      # extract those confluences where every confluence has a LOWER strahler order
      # = first level confluences
      next.df <- m[next.pixel %in% conf,]
      # for each confluence, pick only the confluences that have merging streams of a lower strahler order
      next.df[, max.ord := max(strah_ord), by = .(next.pixel)]
      # first level confluences = merges only strahler orders smaller
      first.level <- unique(next.df[max.ord < s,]$next.pixel)
      # second level confluences = merges strahler orders smaller and same
      sec.level <- unique(next.df[max.ord == s,]$next.pixel)
      
      print(paste0("Starting parallel loop for stream order ", s," and first level streams"))
      temp <-
        foreach(i = 1:length(first.level), .combine = rbind) %dopar% {

          if (m[m$pointid == first.level[i], ]$flag_conf == 1 &
              m[m$pointid == first.level[i], ]$strah_ord == s) {
            # if pointid is a confluence and matches the given stream order
            z <- 1
            k <- first.level[i] # get pointid
            temp.df <- data.frame()
            temp.df[z, "pointid"] <- k
            
            # get pointid of pixels flowing into confluence
            conf.df <- m[m$flag_mouth == 0 & m$next.pixel == k, ]
            conf.df <- conf.df[!is.na(wrt_min.wb),]
              
              # if the confluences are empty (= flooded by lake),
              # then track them down until next confluence
              
              # weigh travel time by discharge for pixels flowing into the confluence
              # and take the sum of all weighted travel time of pixels flowing into confluence
              left.points <-
                cumtime[cumtime$pointid %in% conf.df$pointid, ]$pointid
              # keep only confluences that have travel time
              left.points <-
                left.points[!(is.na(cumtime[cumtime$pointid %in% left.points, ]$travel.time))]
              
              # if WRT is available
              if(length(left.points) > 0L){
                ctime <-
                  sum(cumtime[cumtime$pointid %in% left.points, ][order(pointid),]$travel.time *
                       m[m$pointid %in% left.points, ][order(pointid),]$dis_m3s, na.rm = T)
                # take the sum of discharge, of all pixels flowing into confluence
                cda <- sum(m[m$pointid %in% left.points,]$dis_m3s)
                
                # assign new value to cumulative time in confluence pixel based on:
                # ctime = sum of discharge weighted travel time of pixels flowing into the confluence
                # cda = sum of discharge of pixels flowing into the confluence
                # divided by the residence time within the confluence pixel
                # nx.trav <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
                
                temp.df[z, "travel.time"] <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
      
                # Once we have corrected the confluence travel time by the discharge of the merging streams,
                # we can trail the stream down until the next confluence
                nx.px <- m[m$pointid == k,]$next.pixel
                
                while (m[m$pointid == nx.px, ]$flag_conf != 1) {
                  # while next pixel is not a confluence, do...
                  # calculate cumulative travel time by...
                  if (m[m$pointid == k, ]$flag_mouth != 1) {
                    #print(paste0("I'm at...", k))
                    z <- z + 1
                    # calculate cumulative travel time by...
                    temp.df[z, "pointid"] <- nx.px
                    temp.df[z, "travel.time"] <- 
                      temp.df[temp.df$pointid == k, ]$travel.time + # Time A
                      (m[m$pointid == nx.px, ]$wrt_min.wb * # Time B
                         # adding wrt of next pixel to the weighted wrt of the current pixel
                         (1 - (m[m$pointid == nx.px, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
                            m[m$pointid == nx.px, ]$dis_m3s))
                    
                    # once this is done, overwrite k with the next pixel
                    k <-
                      m[m$pointid == k, ]$next.pixel # so the loop jumps to the next pixel
                    nx.px <- m[m$pointid == k,]$next.pixel
                  } else {
                    #if it hits the mouth, stop
                    break
                  }
                  # if(m[m$pointid == k, "flag_conf"] == 1){
                  #   cat(paste0("\rReached confluence of stream of order ", s,". ", length(conf) - i, " streams within order left."))
                  # }
                  
                }
                
                temp.df <- temp.df[!is.na(temp.df$travel.time), ]
                return(temp.df)
              
                } else {
            # if no stream flowing into confluence has WRT, skip
            break
          } 
        }
      }
      # save results
      cumtime <-
        cumtime[temp, travel.time := i.travel.time, on = .(pointid)]
      
      
      print("Starting while loop to fill in all other secondary level rivers")  
      while(length(sec.level) > 0L){
      # get next.pixels that have travel time already
      prior <- m[m$next.pixel %in% sec.level,]
      prior <- prior[cumtime, travel.time := i.travel.time, on = .(pointid)][!is.na(wrt_min.wb),]
      prior <- prior[flag_main.chan == 1,]
      # calculate how many streams go in each  confluence
      prior[, n.conf := .N, by = .(next.pixel)]
      #any(prior$strah_ord > s)
      # calculate how many of those have already a travel time assigned
      prior[!is.na(travel.time), n.filled := .N, by = .(next.pixel)]
      prior[order(next.pixel),]
      prior <- prior[!is.na(n.filled),][, .(n.conf = unique(n.conf),
                                        n.filled = unique(n.filled)), by = .(next.pixel)]
      prior <- prior[n.filled == n.conf,]$next.pixel
      first <- sec.level[sec.level %in% prior]
      
      print(paste0("Starting parallel loop for stream orders ", s, " to fill in streams of second level"))
      temp <-
        foreach(i = 1:length(first), .combine = rbind) %dopar% {
          #.options.snow = opts
          #for(i in 1:length(first)){
          if (m[m$pointid == first[i], ]$flag_conf == 1 &
              m[m$pointid == first[i], ]$strah_ord == s) {
            # if pointid is a confluence and matches the given stream order
            z <- 1
            k <- first[i] # get pointid
            temp.df <- data.frame()
            temp.df[z, "pointid"] <- k
            
            # get pointid of pixels flowing into confluence
            conf.df <- m[m$flag_mouth == 0 & m$next.pixel == k, ]
            conf.df <- conf.df[!is.na(wrt_min.wb),]
            
            # if the confluences are empty (= flooded by lake),
            # then track them down until next confluence
            
            # weigh travel time by discharge for pixels flowing into the confluence
            # and take the sum of all weighted travel time of pixels flowing into confluence
            left.points <-
              cumtime[cumtime$pointid %in% conf.df$pointid, ]$pointid
            # keep only confluences that have travel time
            left.points <-
              left.points[!(is.na(cumtime[cumtime$pointid %in% left.points, ]$travel.time))]
            
            # if WRT is available
            if(length(left.points) > 0L){
              ctime <-
                sum(cumtime[cumtime$pointid %in% left.points, ][order(pointid),]$travel.time *
                      m[m$pointid %in% left.points, ][order(pointid),]$dis_m3s, na.rm = T)
              # take the sum of discharge, of all pixels flowing into confluence
              cda <- sum(m[m$pointid %in% left.points,]$dis_m3s)
              
              # assign new value to cumulative time in confluence pixel based on:
              # ctime = sum of discharge weighted travel time of pixels flowing into the confluence
              # cda = sum of discharge of pixels flowing into the confluence
              # divided by the residence time within the confluence pixel
              # nx.trav <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
              
              temp.df[z, "travel.time"] <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
            
              # Once we have corrected the confluence travel time by the discharge of the merging streams,
              # we can trail the stream down until the next confluence
              nx.px <- m[m$pointid == k,]$next.pixel
              
              while (m[m$pointid == nx.px, ]$flag_conf != 1) {
                # while next pixel is not a confluence, do...
                  z <- z + 1
                  # calculate cumulative travel time by...
                  temp.df[z, "pointid"] <- nx.px
                  temp.df[z, "travel.time"] <- 
                    temp.df[temp.df$pointid == k, ]$travel.time + # Time A
                    (m[m$pointid == nx.px, ]$wrt_min.wb * # Time B
                       # adding wrt of next pixel to the weighted wrt of the current pixel
                       (1 - (m[m$pointid == nx.px, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
                          m[m$pointid == nx.px, ]$dis_m3s))
              
                  # once this is done, overwrite k with the next pixel
                  k <-
                    m[m$pointid == k, ]$next.pixel # so the loop jumps to the next pixel
                  nx.px <- m[m$pointid == k,]$next.pixel
                  
                  if(is.na(nx.px)){
                    break
                  }
              }
              
              temp.df <- temp.df[!is.na(temp.df$travel.time), ]
              return(temp.df)

            } else {
              # if no stream flowing into confluence has WRT, skip
              break
            }
          }
        }
      # save results
      cumtime <-
        cumtime[temp, travel.time := i.travel.time, on = .(pointid)]
      
      # t <- merge(m, cumtime, all = TRUE)
      # ggplot(t %>% filter(strah_ord <= 2), aes(x = fllength_dec, y = log10(travel.time), colour = strah_ord)) +
      #   geom_point()
      sec.level <- sec.level[!(sec.level %in% first)]
      }
      
      
      print(paste("Finished stream orders...", s))
      # save before moving on to next order
      
    }
    # save final output
    saveRDS(cumtime, paste0("./Objects/cumtime_", y, "-", mon, "_", Sys.Date(), ".rds"))
    rm(cumtime, streamnet, path, m, temp, sources)
  }}
