library(tidyverse)
###########  Relative Abundance plot for top n taxonomy
Rel.plot <- function (rabund,  taxrank, topn, sn, colors ){
  samples_df <- sample_data(rabund)
  samples_df$sample0 = row.names(samples_df)
  samples_df <- as_tibble(samples_df)
  toptaxa = names(head(sort(rowSums(otu_table(rabund)[,sn]), decreasing = TRUE), topn)) 
  rb = prune_samples(sample_names(rabund) %in% sn, rabund)
  rb = prune_taxa(toptaxa, rb)
  df = cbind(data.frame(ID = c(taxa_names(rb), "Other"),
                        taxa = c(tax_table(rb)[, taxrank ], "Other")),
             rbind(otu_table(rb), 1 - colSums(otu_table(rb))))
  
  df = df[match(c("Other",toptaxa ), df$ID),]
  df = df %>% mutate_at(vars(taxa), factor) %>% mutate(taxa = factor(taxa, levels = df$taxa))
  prel <- pivot_longer(df, cols = -c(ID, taxa), names_to = "sample0", values_to = "Abundance") %>%
    left_join(samples_df %>% select(sample ,type3, sample0)) %>%
    
    ggplot(aes(sample,  Abundance, fill = taxa)) + geom_col(width=0.7, color = "grey") +
    scale_fill_manual(values = colors)+ scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    labs( y = "Relative abundance", x= " ")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text.x = element_text(size = 12),
          strip.background =element_rect(fill="white")) +
    facet_wrap(~type3, 1, scales = "free_x")+
    guides(fill=guide_legend(title= taxrank ))
  prel
}

########### beta diversity plot
betadiv_plot <- function(phydat, nrow, scale, title ) {
  immc.ord <- ordinate(phydat, "NMDS", "bray")
  p <- plot_ordination(phydat, immc.ord,  color="group", title = title ) +  
    geom_point( size=3)+
    theme_bw() 
  
  p <- p  + scale_color_manual("Method ", values =c("darkcyan", "darkgoldenrod2", "#FC4E07"))+
    scale_shape_manual(" ", values=c(16,  17, 15), labels = c("All", "+control" ))+
    geom_text(mapping = aes(label = F), size = 3, hjust = -0.1, vjust=1.5) +
    theme(strip.text.x = element_text(size = 12),
          strip.background =element_rect(fill="white"))+
    facet_wrap(~type2, nrow, scales= scale ) 
  p
}

########### Abundance/Relative Abundance Plot

abund_plot <- function(abunddat, x.var, y.var, fill.var, x=0, y, col ){
  
  abunddat$y.var <- abunddat[, y.var]
  
  p1=abunddat %>% filter (environment=="Lab_Q")%>%
    ggplot( aes_string(x=x.var, y=y.var, fill=fill.var)) +
    geom_bar(stat="identity", position="stack", width=0.5)+
    scale_fill_manual("",values = col)+
    theme_minimal() + xlab("")+ ylab("") + 
    guides(size = guide_legend(title ="Archaea"))+
    theme(legend.position = "none",
          axis.text.x=element_blank())+
    facet_grid(~factor(method, levels=c("Coral: Qiagen (Lab)", "Mangrove: Qiagen (Lab)")), scales = "free_x") 
  
  if (sum(abunddat$y.var) > nrow(abunddat)) {
    p1 <- p1+
      scale_y_continuous(limits =c(x, y)  , labels = unit_format(unit = "M", scale = 1e-6))
  }else {
    p1 <- p1+
      scale_y_continuous(limits =c(x, y))
  }
  
  p2= abunddat %>% filter (environment=="Field") %>%
    ggplot(aes_string(x=x.var, y=y.var,  fill=fill.var)) +
    geom_bar(stat="identity", position="stack", width=0.5)+
    scale_fill_manual("",values = col)+
    theme_minimal() + xlab("")+ ylab("Abundance") + 
    guides(size = guide_legend(title ="Archaea"))+
    theme(legend.position = "none",
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.x=element_blank())+
    facet_grid(~factor(method, levels = c("xTitan (Boat)" , "xTitan (In situ)")), scales = "free_x")
  
  if (sum(abunddat$y.var) > nrow(abunddat)) {
    p2 <- p2+
      scale_y_continuous(limits =c(x, y)  , labels = unit_format(unit = "M", scale = 1e-6))
  }else {
    p2 <- p2+
      scale_y_continuous(limits =c(x, y))
  }
  
  p3=abunddat %>% filter (environment=="Lab") %>%
    ggplot(aes_string(x=x.var, y=y.var, fill=fill.var)) +
    geom_bar(stat="identity", position="stack", width=0.5)+
    scale_fill_manual("",values = col)+
    theme_minimal() + xlab("")+ ylab("")+
    guides(size = guide_legend(title ="Archaea"))+
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_grid(~factor(method, levels = c("xTitan (Lab)" , "X-Titan (Lab)"), 
                       labels= c("xTitan (Lab)" , "xTitan (Lab) ")), scales = "free_x")
  
  if (sum(abunddat$y.var) > nrow(abunddat)) {
    p3 <- p3+
      scale_y_continuous(limits =c(x, y)  , labels = unit_format(unit = "M", scale = 1e-6))
  }else {
    p3 <- p3+
      scale_y_continuous(limits =c(x, y))
  }
  
  gridExtra::grid.arrange(p1,p2,p3, heights=c(1,1,1.4))
  
}

################# Chord diagram
chord <- function(mat, grid.col,orders, col, label = TRUE) {
  rn = rownames(mat)
  cn = colnames(mat)
  chordDiagram(mat, 
               order = orders,
               grid.col=grid.col, 
               annotationTrack = c("grid", "axis"),
               row.col =col,
               big.gap = 15, small.gap = 2,
               transparency = 0.3,
               link.zindex = rank(mat),
               preAllocateTracks = list(track.height = 0.3), annotationTrackHeight = uh(3, "mm"),
               diffHeight = mm_h(3),
               directional = 1, direction.type = "diffHeight")
  
  if (label == "TRUE") {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      label = data.table::tstrsplit(sector.name, "_")[[1]]
      if(sector.name %in% rn) {
        circos.text(mean(xlim), ylim[1], label, facing = "reverse.clockwise",  col="blue",adj = c(1, 1.5),cex=0.9)
      }
      if(sector.name %in% cn) {
        circos.text(mean(xlim), ylim[1]+0.2, label, facing = "reverse.clockwise", adj = c(0.7, 1), col="red",cex=0.9)
      }
    }, bg.border = NA)
  }
  circos.clear()
}

