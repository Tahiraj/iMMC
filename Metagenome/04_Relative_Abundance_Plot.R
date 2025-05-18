library(phyloseq)
library(RColorBrewer)
library(wesanderson)

rm(list = ls())
immc_decontam <- readRDS("immc_decontam.rds")
source("functions.R")
#==================================================================
# Relative abundace figures 
#==================================================================

Rel.plot0 <- function (rabund, taxrank, topn, sn ){
  
  sn = sample_names(sample_data(rabund)[sample_data(rabund)$type2=="Coral1" |sample_data(rabund)$type2=="Coral2" ,])
  pg_coral <- Rel.plot (rabund,  taxrank , topn, sn, colors=colors_coral )
  
  sn = sample_names(sample_data(rabund)[sample_data(rabund)$type02=="Coral1-W"|sample_data(rabund)$type02=="Coral2-W" ,])
  pg_coralw <- Rel.plot (rabund,  taxrank , topn, sn, colors=colors_coralw)
  
  sn = sample_names(sample_data(rabund)[sample_data(rabund)$type== "Mangrove: Water",])
  pg_water <- Rel.plot (rabund,  taxrank , topn, sn, colors=colors_water )
  
  sn = sample_names(sample_data(rabund)[sample_data(rabund)$type== "Mangrove: Sediment",])
  pg_sediment <- Rel.plot (rabund,  taxrank, topn, sn, colors=colors_sediment )
  
  gridExtra::grid.arrange(pg_coral, pg_coralw, pg_water, pg_sediment, nrow=2, widths=c(3,  0.5, 0.5, 2.2), 
                          layout_matrix = rbind(c(1,1,NA,2), c(3, NA,4,4)))
}

colors_coral <- c("gray","#E41A1C", "#596A98", "#449B75", "#6B886D", "#AC5782", "#FF7F00", "#FFE528", "#C9992C", "#C66764", "#E485B7", "black") 
colors_coralw <- c("gray", rev(colorRampPalette(brewer.pal(8,"Spectral")) (10)))
colors_water <- c("gray", colorRampPalette(brewer.pal(10,"PiYG")) (10))
colors_sediment <- c("gray", wes_palette("Cavalcanti1", 10, type = "continuous"))

gp_g = tax_glom(immc_decontam, taxrank = "Genus") 
rabundg = transform_sample_counts(gp_g, function(x) {x/sum(x)} )
Rel.plot0 (rabundg, taxrank = "Genus", topn = 10, sn = sn )

#===============================================================================
gp_p = tax_glom(immc_decontam, taxrank = "Phylum")
rabundp = transform_sample_counts(gp_p, function(x) {x/sum(x)} )
Rel.plot0 (rabund= rabundp, taxrank = "Phylum", topn = 10, sn = sn )

#===============================================================================
gp_p = tax_glom(immc_decontam, taxrank = "Family")
rabundf = transform_sample_counts(gp_p, function(x) {x/sum(x)} )
Rel.plot0 (rabund= rabundf, taxrank = "Family", topn = 10, sn = sn )
#===============================================================================

pdf("Relative_Abundance_decontam.pdf", width=18, height=14)
Rel.plot0 (rabund=  rabundg, taxrank = "Genus", topn = 10, sn = sn )
Rel.plot0 (rabund= rabundf, taxrank = "Family", topn = 10, sn = sn )
Rel.plot0 (rabund= rabundp, taxrank = "Phylum", topn = 10, sn = sn )
dev.off()


