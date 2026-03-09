# Collection of functions for smooth plotting ----------------------------------------


theme_custom <- function(base_size = 11, base_family = "", base_line_size = base_size/22, 
                         base_rect_size = base_size/22, half_line = base_size/2){
  
  theme_grey(base_size = base_size, base_family = base_family,
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA),
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = 'grey20', colour = "grey20"),
          strip.text = element_text(colour = 'white',
                                    size = rel(0.8), margin = margin(0.8 * half_line, 
                                                                     0.8 * half_line, 0.8 * half_line, 0.8 * half_line)),
          legend.key = element_rect(fill = "white", colour = NA), complete = T)
}

# Custom theme with smaller legends
theme_cust <- function(base_theme = "bw", base_size = 11, half_line = base_size/2,
                       border = F){
  if(base_theme == "bw"){
    t <- theme_bw(base_size = base_size)
  }
  
  if(base_theme == "pubr"){
    if(border == T){
      t <- theme_pubr(base_size = base_size, border = T)
    } else {
      t <- theme_pubr(base_size = base_size)
    }
  }
  
  t %+replace% theme(legend.position = "right", 
                     legend.spacing = unit(half_line / 2, "pt"), # spacing between legends
                     legend.key.size = unit(0.7, "lines"), # size of legend symbol box
                     legend.box.spacing = unit(1.5 * half_line, "pt"), # spacing between legend and plot
                     legend.text = element_text(size = unit(base_size - 4, "pt")), 
                     legend.title = element_text(size = unit(base_size - 3, "pt"), hjust = 0)
  )
  
}

plot_pca<- function(pca, meta, by, axes = 1:2, scale){
  # by = "OTU" or "MF"
  
  loadings  <- pca$rotation
  scores <- pca$x
  var <- summary(pca)$importance
  
  loadings <- setDT(data.frame(loadings), keep.rownames = by)
  setDT(meta)
  
  loadings <- loadings[meta, , on = c(by)]
  loadings[ns.s == "n.s.", sAU := "n.s."]
  loadings[, sAU := factor(sAU, levels = c("increase", "non-linear increase",
                                           "unimodal",
                                           "non-linear decrease", "decrease","n.s."),
                           labels = c("Increase", "Non-linear increase",
                                      "Unimodal",
                                      "Non-linear decrease", "Decrease", "n.s."))]
  
  scores <- data.frame(scores)
  scores$sites <- row.names(scores)
  scores <- separate(scores, sites, into = c("Year","Season", "Site"), sep = "[.]")
  scores$Site <- as.numeric(scores$Site)
  scores$next.x <- lead(scores[, axes[1]], 1)
  scores$next.y <- lead(scores[, axes[2]], 1)
  
  cols <- c("#B2182B", "#F4A582",
            "#F7F7F7", "#92C5DE", "#2166AC", "black")
  names(cols) <-  c("Increase", "Non-linear increase",
                    "Unimodal",
                    "Non-linear decrease", "Decrease", "n.s.")
  
  # colour arrows and sites 
  # cols <- c('#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', 
  #                    '#F67E4B', '#DD3D2D') # super red '#A50026', bad vals:'#FFFFFF'
  #                    scale_colour_gradientn(colours = cols, trans = "reverse")
  ggplot()+
    theme_bw() +
    #coord_fixed() +
    geom_point(data = loadings[sAU == "n.s.",], aes(x = .data[[paste0("PC", axes[1])]]*scale,
                                                    y = .data[[paste0("PC", axes[2])]]*scale),
               colour = "grey70",  size = 1.5) +  
    geom_point(data = loadings[ns.s == "sign",], aes(x = .data[[paste0("PC", axes[1])]]*scale, 
                                                     y = .data[[paste0("PC", axes[2])]]*scale, fill = sAU),
               size = 2.5, shape = 21) + #, label = OTU
    geom_segment(data = scores, aes(x = (.data[[paste0("PC", axes[1])]] + next.x)/2, xend = next.x,
                                    y = (.data[[paste0("PC", axes[2])]] + next.y)/2, yend = next.y)) +
    geom_segment(data = scores, aes(x = .data[[paste0("PC", axes[1])]], 
                                    xend = (.data[[paste0("PC", axes[1])]] + next.x)/2, 
                                    y = .data[[paste0("PC", axes[2])]], 
                                    yend = (.data[[paste0("PC", axes[2])]] + next.y)/2),
                 arrow = arrow(length = unit(0.3, "cm")), show.legend=FALSE) +
    geom_point(data = scores, aes(x = .data[[paste0("PC", axes[1])]], 
                                  y = .data[[paste0("PC", axes[2])]]), size = 3) +
    labs(x = paste0("PC",  axes[1]," (", round(var[2, axes[1]] * 100,
                                               digits = 1),"%)"),
         y = paste0("PC",  axes[2]," (", round(var[2, axes[2]] * 100,
                                               digits = 1),"%)"),
         title = paste(unique(loadings$dataset),":", unique(loadings$Year), unique(loadings$Season), sep = " ")) +
    scale_fill_manual(values = cols[names(cols) %in% unique(loadings$sAU)], name = "Spatial pattern") 
}

plot_pcoa <- function(pcoa, axes = 1:2){
  pdataframe <- data.frame(Sample = as.character(row.names(pcoa$vectors)),
                           pcoa$vectors[,axes],
                           stringsAsFactors = F)  # extract site scores
  pb.var <- data.frame(Axes = axes,
                       var = round(100 * pcoa$values$Eigenvalues[axes] / sum(pcoa$values$Eigenvalues), 2),
                       stringsAsFactors = F)
  
  #pdataframe <- merge(pdataframe, meta, by = "Sample")
  #pdataframe$Sample <- as.character(pdataframe$Sample)
  
  (p <- ggplot(pdataframe, aes_string(x = paste0("Axis.", axes[1]),
                                      y = paste0("Axis.", axes[2]))) +
      theme_bw() +
      geom_hline(yintercept =  0, colour = "grey80", size = 0.4) +
      geom_vline(xintercept = 0, colour = "grey80", size = 0.4) +
      geom_point() +
      coord_fixed(1) + # ensure aspect ratio
      labs(x = paste0("PC",  axes[1]," (", round(pb.var[pb.var$Axes == axes[1],"var"],
                                                 digits = 1),"%)"),
           y = paste0("PC",  axes[2]," (", round(pb.var[pb.var$Axes == axes[2],"var"],
                                                 digits = 1),"%)")) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  )
}

# originals from theme_grey
#legend.key.size = unit(1.5, "lines") # size of legend symbol box
#legend.spacing = unit(2 * half_line, "pt") # spacing between legends
#legend.box.spacing = unit(2 * half_line, "pt") # spacing between legend and plot
#plot.title = element_text(size = rel(1.2)