library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library("metagMisc") # Convert phyloseq object to data frame 

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
#pdf("Relative_Abundance_decontam.pdf", width=18, height=14)
Rel.plot0 (rabund=  rabundg, taxrank = "Genus", topn = 10, sn = sn )
Rel.plot0 (rabund= rabundf, taxrank = "Family", topn = 10, sn = sn )
Rel.plot0 (rabund= rabundp, taxrank = "Phylum", topn = 10, sn = sn )
#dev.off()

#================================================================================
### Coral Abundance plot for species within selected taxonomy (Family)
#==================================================================================
samples_df0 <- data.frame(sample_data(immc_decontam))  
samples_df0$sample <- rownames(samples_df0)

family <- c("Endozoicomonadaceae", "Spiroplasmataceae", "Fulvivirgaceae", "Francisellaceae")

rb_sp = transform_sample_counts(immc_decontam, function(x) {x/sum(x)} )

sn = sample_names(sample_data(rb_sp)[sample_data(rb_sp)$type2=="Coral1"|sample_data(rb_sp)$type2=="Coral2" , ])

rb0 = prune_samples(sample_names(rb_sp) %in% sn, rb_sp)

#===========================================================
rb =  subset_taxa(rb0, Family=="Endozoicomonadaceae")

topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 6)) 

rb = prune_taxa(topn, rb)
rb
df1 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df1 = df1[topn,]
df1 = df1 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df1$Species))

#=================================================================================

rb =  subset_taxa(rb0, Family=="Spiroplasmataceae")

topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)
rb
df2 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df2 = df2[topn,]
df2 = df2 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df2$Species))

#=================================================================================

rb =  subset_taxa(rb0, Family=="Fulvivirgaceae")

topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)

df3 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df3 = df3[topn,]

df3 = df3 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df3$Species))

#================================================================================

rb =  subset_taxa(rb0, Family=="Francisellaceae")

topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)
rb
df4 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df4 = df4[topn,]

df4 = df4 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df4$Species))

df <- bind_rows(df1, df2, df3, df4)

df %>% group_by(Family ) %>% tally()

colors_coral <-c(
  colorRampPalette(c( "darkgreen", "olivedrab1", "whitesmoke")) (6),
  colorRampPalette(c("#596A98",  "lightblue", "azure", "white")) (5),
  colorRampPalette(c( "orange", "yellow","lightyellow")) (5),
  colorRampPalette(c( "#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF")) (5))

differential_Abundance_Coral <- df %>% 
  pivot_longer( cols = -c(ID, Species, Family ), names_to = "sample", values_to = "Abundance") %>%
  left_join(samples_df0 %>% select(sample ,type3)) %>% 
  filter(type3 %in% c("Coral1","Coral2" )) %>%
  
  ggplot(aes(sample,  Abundance, fill = Species )) + geom_col(width=0.7, color = "grey") +
  scale_fill_manual(values = colors_coral)+ scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs( y = "", x= " ")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background =element_rect(fill="white")) +
  guides(fill = guide_legend(ncol = 1))+
  facet_grid(Family ~ type3, scales = "free_x")


pdf("differential_Abundance_Coral.pdf", width=7, height=10)
differential_Abundance_Coral
dev.off()
#================================================================================
#Mangrove
#================================================================================
sn = sample_names(sample_data(rb_sp)[sample_data(rb_sp)$type2=="Mangrove: Water"|sample_data(rb_sp)$type2=="Mangrove: Sediment" , ])

rb0 = prune_samples(sample_names(rb_sp) %in% sn, rb_sp)

family=c("Balneolaceae", "Pirellulaceae", "Flavobacteriaceae", "Rhodobacteraceae")

#================================================================================

rb =  subset_taxa(rb0, Family=="Balneolaceae")
rb
topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)

df1 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df1 = df1[topn,]
df1 = df1 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df1$Species))

#=================================================================================

rb =  subset_taxa(rb0, Family=="Pirellulaceae")
rb
topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)

df2 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df2 = df2[topn,]
df2 = df2 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df2$Species))

#=================================================================================

rb =  subset_taxa(rb0, Family=="Flavobacteriaceae")
rb
topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)

df3 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df3 = df3[topn,]

df3 = df3 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df3$Species))

#================================================================================

rb =  subset_taxa(rb0, Family=="Rhodobacteraceae")
rb
topn = names(head(sort(rowSums(otu_table(rb)[,sn]), decreasing = TRUE), 5)) 

rb = prune_taxa(topn, rb)

df4 = cbind(data.frame(ID = taxa_names(rb)),
            tax_table(rb)[,c("Family","Species")],
            otu_table(rb))

df4 = df4[topn,]

df4 = df4 %>% mutate_at(vars(Species), factor) %>% mutate(Species = factor(Species, levels = df4$Species))

df <- bind_rows(df1, df2, df3, df4)

df %>% group_by(Family ) %>% tally()

colors_mangrove <-c(
  colorRampPalette(c( "darkgreen", "olivedrab1", "whitesmoke")) (5),
  colorRampPalette(c("#596A98",  "lightblue", "azure", "white")) (5),
  colorRampPalette(c("orange4","orange2","gold1",  "yellow2","lightyellow")) (5),
  colorRampPalette(c( "#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF")) (5))

differential_Abundance_Mangrove <- df %>% 
  pivot_longer(cols = -c(ID, Species, Family ), names_to = "sample", values_to = "Abundance") %>%
  left_join(samples_df0 %>% filter(environment0 == "Mangrove") %>% select(sample ,type3)) %>%
  
  ggplot(aes(sample,  Abundance, fill = Species )) + geom_col(width=0.7, color = "grey") +
  scale_fill_manual(values = colors_mangrove)+ scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  labs( y = "", x= " ")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background =element_rect(fill="white")) +
  guides(fill = guide_legend(ncol = 1))+
  facet_grid(Family ~ type3, scales = "free_x")

pdf("differential_Abundance_Mangrove.pdf", width=9, height=11)
differential_Abundance_Mangrove
dev.off()

#===================================
#Differential expression analysis
#===================================
#BiocManager::install("DESeq2")
library(DESeq2)
deseq <- function(phy_dat, sn, alpha) {
  rb = prune_samples(sample_names(phy_dat) %in% sn, phy_dat)
  
  dds <- phyloseq_to_deseq2(rb, ~group)
  dds <- DESeq(dds, test="Wald", fitType="parametric")
  
  res = results(dds, cooksCutoff = FALSE)
  
  alpha = alpha 
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy_dat)[rownames(sigtab), ], "matrix"))
  return (sigtab )
}

sigtab_taxa <- function(sigtab, phy_dat, sn,a,b ){
  taxa <- rownames(sigtab )
  rb_M = prune_samples(sample_names(phy_dat) %in% sn, phy_dat)
  rb_M = prune_taxa(taxa_names(rb_M) %in% taxa, rb_M)
  rb_df <- phyloseq_to_df(rb_M, addtot=T)
  rb_df0 <- rb_df %>%
    rowwise() %>% 
    mutate(Score1 = sum(c_across(starts_with(a))), Score2 = sum(c_across(starts_with(b)))) %>% 
    select(OTU, Score1,Score2, Total ) %>%
    tibble::column_to_rownames("OTU")
  
  sigtab <- merge(sigtab, rb_df0,  by = 'row.names', all = TRUE) %>%
    mutate(sp=str_remove(`Row.names`, "sp" )) %>%
    arrange(as.numeric(sp)) %>% select(-sp)
  
  return(sigtab)
}


deseqtab <- function (phy_dat) {
  
  #===================================================================================
  # Mangrov-Water: 
  #===================================================================================
  # xTitan(L) vs xTitan(F)
  # phy_dat = gp_p
  sn_Mang = sample_names(sample_data(phy_dat)[sample_data(phy_dat)$type2=="Mangrove: Water" ,]) 
  sn = sn_Mang[grep ("titan", sn_Mang)]
  
  sigtab_MW_tt <- deseq (phy_dat, sn, alpha=0.05)
  
  sigtab_MW_tt <- sigtab_MW_tt %>% mutate (Biome= "Mangrove: Water", Method= "xTitan(L):xTitan(F)")
  
  sigtab_MW_tt0 <- sigtab_taxa (sigtab=sigtab_MW_tt, phy_dat , sn, a="titanL", b="titanF" )
  
  #====================================
  #Qiagen Vs  xTitan(L) 
  
  sn = sn_Mang[grep ("titanL|qiagen", sn_Mang)]
  
  sigtab_MW_tq <- deseq (phy_dat, sn, alpha=0.05)
  sigtab_MW_tq <- sigtab_MW_tq %>% mutate (Biome= "Mangrove: Water", Method= "xTitan(L):Qiagen") 
  
  sigtab_MW_tq0 <- sigtab_taxa (sigtab=sigtab_MW_tq, phy_dat, sn, a="titanL", b="Qiagen" )
  
  #===================================================================================
  # Mangrov-Sediment: 
  #===================================================================================
  # xTitan(L) vs xTitan(F)
  
  sn_Mang = sample_names(sample_data(phy_dat)[sample_data(phy_dat)$type2=="Mangrove: Sediment" ,]) 
  sn = sn_Mang[grep ("titan", sn_Mang)]
  
  sigtab_MS_tt <- deseq (phy_dat, sn, alpha=0.05)
  sigtab_MS_tt <- sigtab_MS_tt %>% mutate (Biome= "Mangrove: Sediment", Method= "xTitan(L):xTitan(F)")
  
  sigtab_MS_tt0 <- sigtab_taxa (sigtab=sigtab_MS_tt, phy_dat, sn, a="titanL", b="titanF" )
  #================================================================
  #Qiagen Vs  xTitan(L) 
  
  sn = sn_Mang[grep ("titanL|qiagen", sn_Mang)]
  
  sigtab_MS_tq <- deseq (phy_dat, sn, alpha=0.05)
  sigtab_MS_tq <- sigtab_MS_tq %>% mutate (Biome= "Mangrove: Sediment", Method= "xTitan(L):Qiagen")
  
  sigtab_MS_tq0 <- sigtab_taxa (sigtab=sigtab_MS_tq, phy_dat, sn, a="titanL", b="Qiagen" )
  
  sigtabM <-bind_rows(sigtab_MW_tt0, sigtab_MW_tq0, sigtab_MS_tt0, sigtab_MS_tq0)
  
  return (sigtabM)  
}

 #================================================================
sigtabM <- deseqtab (immc_decontam )

#write.csv (sigtabM, "diffabund_ASV_sigtabM.csv", row.names = FALSE)
dim(sigtabM)
sigtabw <- sigtabM %>% 
  mutate(Group=ifelse(log2FoldChange > 0,"A","B")) %>%
  filter(log2FoldChange > 10 | log2FoldChange < -10) %>%
  filter(Biome == "Mangrove: Water")

p1 <- sigtabw %>% filter(Method =="xTitan(L):xTitan(F)") %>%
  ggplot(aes(x = reorder(Species, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0))+
  geom_bar(stat = "identity", width=0.5, position=position_dodge2(preserve = "single"), aes(fill=Group) )+
  scale_fill_viridis_d( option = "D", begin = 0.2, end = 0.8, direction = 1) + 
  labs (title="Mangrove: Water" , x ="Species")+
  guides(fill=guide_legend(title=""))+
  coord_flip()
p1 <- p1  +theme_bw()+theme(legend.position = "none")
p1

p11 <- sigtabw %>% filter(Method =="xTitan(L):Qiagen") %>%
  ggplot(aes(x = reorder(Species, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0))+
  geom_bar(stat = "identity", width=0.5, position=position_dodge2(preserve = "single"), aes(fill=Group) )+
  scale_fill_viridis_d( option = "D", begin = 0.2, end = 0.8, direction = 1) + 
  labs (title="Mangrove: Water" , x ="")+
  guides(fill=guide_legend(title=""))+
  coord_flip()
p11 <- p11+ theme_bw()+ theme(legend.position = "none")
p11

sigtabm <- sigtabM %>% 
  mutate(Group=ifelse(log2FoldChange > 0,"A","B")) %>%
  filter(log2FoldChange > 10 | log2FoldChange < -10) %>%
  #  filter(Method=="xTitan(L):xTitan(F)") %>%
  filter(Biome == "Mangrove: Sediment")

p2 <- sigtabm %>% filter(Method =="xTitan(L):xTitan(F)") %>%
  ggplot(aes(x = reorder(Species, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0))+
  geom_bar(stat = "identity",width=0.5,  position=position_dodge2(preserve = "single"), aes(fill=Group) )+
  scale_fill_viridis_d( option = "D", begin = 0.2, end = 0.8, direction = 1) +
  labs (title="Mangrove: Sediment" , x="Species")+
  guides(fill=guide_legend(title=""))+
  coord_flip()
p2 <- p2 + theme_bw()+ theme(legend.position = "none")

p22 <- sigtabm %>% filter(Method =="xTitan(L):Qiagen") %>%
  ggplot(aes(x = reorder(Species, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0))+
  geom_bar(stat = "identity",width=0.5,  position=position_dodge2(preserve = "single"), aes(fill=Group) )+
  scale_fill_viridis_d( option = "D", begin = 0.2, end = 0.8, direction = 1) +
  labs (title="Mangrove: Sediment" , x="")+
  guides(fill=guide_legend(title=""))+
  coord_flip()
p22 <- p22 + theme_bw()+ theme(legend.position = "none")

grid.arrange(p2 ,p22, ncol=2, heights = c(1.5, 1), layout_matrix = rbind(c(1,2), c(1, NA)))
grid.arrange(p1 ,p11, ncol=2, heights = c(1, 1,1), layout_matrix = rbind(c(1,2), c(1, NA),  c(NA, NA)))

pdf("Deseq_log2Fold.pdf", width=10, height=10)
grid.arrange(p2 ,p22, ncol=2, heights = c(1.5, 1), layout_matrix = rbind(c(1,2), c(1, NA)))
grid.arrange(p1 ,p11, ncol=2, heights = c(1, 1,1), layout_matrix = rbind(c(1,2), c(1, NA),  c(NA, NA)))

dev.off()
#=================================================================
abund <- read_csv("diffabund_ASV_sigtabM.csv")

tab <- abund %>%  filter(grepl('Mang', Biome )) %>%
  group_by(Biome, Method) %>%
  tally()

abund %>%  filter(grepl('Mang', Biome )) %>%
  select(Biome, Method) %>%
  table() %>%
  chisq.test() 
# Biome and Methods are correlated

abund %>%  filter(grepl('Mang', Biome )) %>%
  select(Biome, Method) %>%
  ggplot() +
  aes(x = Biome, fill = Method) +
  geom_bar(position = "dodge", width=0.7) +xlab("")+ ylab("") +
  theme_bw()

#=================================================================


