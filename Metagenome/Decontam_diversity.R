library(readxl) 
library(pairwiseAdonis)

rm(list = ls())

immc_decontam <- readRDS("immc_decontam.rds")
source("functions.R")

#==================================================================
# Diversity plots
#================================================= Beta-diversity
#====================== Mangrove

sn = sample_names(sample_data(immc_decontam)[sample_data(immc_decontam)$type0=="Mangrove",])
immc_mangrov = prune_samples(sample_names(immc_decontam) %in% sn, immc_decontam)
div_mang <- betadiv_plot (phydat=immc_mangrov, nrow= 1, scale = "free", title = "Taxonomy")

#====================== Corals

sn = sample_names(sample_data(immc_decontam)[sample_data(immc_decontam)$type0=="Coral",])
immc_coral = prune_samples(sample_names(immc_decontam) %in% sn, immc_decontam)
div_coral <- betadiv_plot (phydat=immc_coral,nrow= 1, scale = "free", title = "Taxonomy")

#pdf("Diversity_plot_taxonomy.pdf", width=12, height=8)
gridExtra::grid.arrange(div_mang, div_coral, widths=c(1,1,0.7),  layout_matrix = rbind(c(1, 1, NA), c(2,2,2)))
#dev.off()
#================================================================================
# Permanova test
#================================================================================
 
OTU_decontaminated <- otu_table(immc_decontam)
OTUp <- t(OTU_decontaminated) 

samples_df <- data.frame(sample_data(immc_decontam))  
samples_dfp <- samples_df %>% select(type2, sample, group )

OTUp <- merge(samples_dfp , OTUp, by = 'row.names', all = TRUE) 

OTUp_Coral1 <- OTUp %>% filter (type2 == "Coral1") %>% select( -c(Row.names, type2, sample, group ))
samples_Coral1 <- samples_dfp %>% filter (type2 == "Coral1")

permanova <- adonis2(OTUp_Coral1~ group ,
                     data = samples_Coral1 , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_Coral1~ group ,
                     data = samples_Coral1 , permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++

OTUp_Coral2 <- OTUp %>% filter (type2 == "Coral2") %>% select( -c(Row.names, type2, sample, group ))
samples_Coral2 <- samples_dfp %>% filter (type2 == "Coral2")

permanova <- adonis2(OTUp_Coral2~ group ,
                     data = samples_Coral2 , permutations=999, method = "bray")

permanova

permanova_p <- pairwise.adonis2(OTUp_Coral2~ group ,
                                data = samples_Coral2 , permutations=999, method = "bray")
permanova_p
#++++++++++++++++++++++++++++++

OTUp_CoralW <- OTUp %>% filter (type2 == "Coral-W") %>% select( -c(Row.names, type2, sample, group ))
samples_CoralW <- samples_dfp %>% filter (type2 == "Coral-W")

permanova <- adonis2(OTUp_CoralW~ group ,
                     data = samples_CoralW , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_CoralW~ group ,
                                data = samples_CoralW, permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++

OTUp_MangroveW <- OTUp %>% filter (type2 == "Mangrove: Water") %>% select( -c(Row.names, type2, sample, group ))
samples_MangroveW  <- samples_dfp %>% filter (type2 == "Mangrove: Water")

permanova <- adonis2(OTUp_MangroveW ~ group ,
                     data = samples_MangroveW  , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_MangroveW~ group ,
                                data = samples_MangroveW, permutations=999, method = "bray")
permanova_p
#++++++++++++++++++++++++++++++

OTUp_MangroveS <- OTUp %>% filter (type2 == "Mangrove: Sediment") %>% select( -c(Row.names, type2, sample, group ))
samples_MangroveS <- samples_dfp %>% filter (type2 == "Mangrove: Sediment")

permanova <- adonis2(OTUp_MangroveS~ group ,
                     data = samples_MangroveS , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_MangroveS~ group ,
                                data = samples_MangroveS, permutations=999, method = "bray")
permanova_p


biome <- c("Coral1", "Coral2", "Coral-W", "Mangrove: Water",  "Mangrove: Sediment")
method <- c("Qiagen (Lab)", "xTitan (Field)", "xTitan (Lab)" )

#=================================================================
# https://chiliubio.github.io/microeco_tutorial/model-based-class.html
# https://www.shermstats.com/2019/02/04/feature-selection-with-boruta-part-2/
  
#======================================================================
# Functional Analysis

#oldnames <- c(C1, C2, C3, C4, C5, C6, C7, C8, C9, M1, M2, M3, M4, M5, M6, M7, M8, M9, PC) 
#newnames <- c(C01, C02, C03, C04, C05, C06, C07, C08, C09, M01, M02, M03, M04, M05, M06, M07, M08, M09, PC1)
abund <- read_excel("KO.Abund.xlsx", sheet ="abundance")
samples_df <- read_excel("KO.Abund.xlsx", sheet ="Sample")%>% 
  filter(!sample%in%c("qiagen_PC1", "titanF_PC3","titanL_PC2")) 
oldnames <- samples_df$F
newnames <- samples_df$sample
data.table::setnames(abund,oldnames,newnames, skip_absent=TRUE)

otu_dat <- abund  %>% 
  select(-PC) %>%
  rowwise() %>% 
  filter(sum(c_across("titanF_C1":"titanF_M9")) != 0) 

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample") 
colnames(samples_df)[1] <- "sample"

tax_dat <- read_excel("KO.Abund.xlsx", sheet ="taxa")

OTU = otu_table(as.matrix(otu_dat[,-1]), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_dat))
samples = sample_data(samples_df)


ko_ps <- phyloseq(OTU, TAX, samples)
taxa_names(ko_ps ) <- otu_dat$kegg
ko_ps <- prune_taxa(!taxa_names(ko_ps )=="Unclassified", ko_ps )


#====================== Mangrove

sn = sample_names(sample_data(ko_ps)[sample_data(ko_ps)$type0 =="Mangrove",])
immc_mangrov = prune_samples(sample_names(ko_ps) %in% sn, ko_ps)
fun_mang <- betadiv_plot (phydat=immc_mangrov, nrow= 1, scale = "free", title = "Kegg: Functional" )

fun_mang
#====================== Corals

sn = sample_names(sample_data(ko_ps)[sample_data(ko_ps)$type0 =="Coral",])
immc_coral =  prune_samples(sample_names(ko_ps) %in% sn, ko_ps)
fun_coral <- betadiv_plot (phydat=immc_coral,nrow= 1, scale = "free", title = "Kegg: Functional")

#pdf("Diversity_plot_functional.pdf", width=12, height=8)
gridExtra::grid.arrange(fun_mang, fun_coral, widths=c(1,1,0.7),  layout_matrix = rbind(c(1, 1, NA), c(2,2,2)))
#dev.off()

#otu_table(ko_ps)
OTUp <- t(otu_table(ko_ps)) 

samples_dfp <- samples_df %>% select(type2, sample, group )

OTUp <- merge(samples_dfp , OTUp, by = 'row.names', all = TRUE) 

#+++++++++++++++++++++++++++++++++==++=====+++++=

OTUp_Coral1 <- OTUp %>% filter (type2 == "Coral1") %>% select( -c(Row.names, type2, sample, group ))
samples_Coral1 <- samples_dfp %>% filter (type2 == "Coral1")

permanova <- adonis2(OTUp_Coral1~ group ,
                     data = samples_Coral1 , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_Coral1~ group ,
                                data = samples_Coral1 , permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++

OTUp_Coral2 <- OTUp %>% filter (type2 == "Coral2") %>% select( -c(Row.names, type2, sample, group ))
samples_Coral2 <- samples_dfp %>% filter (type2 == "Coral2")

permanova <- adonis2(OTUp_Coral2~ group ,
                     data = samples_Coral2 , permutations=999, method = "bray")

permanova

permanova_p <- pairwise.adonis2(OTUp_Coral2~ group ,
                                data = samples_Coral2 , permutations=999, method = "bray")
permanova_p
#++++++++++++++++++++++++++++++

OTUp_CoralW <- OTUp %>% filter (type2 == "Coral-W") %>% select( -c(Row.names, type2, sample, group ))
samples_CoralW <- samples_dfp %>% filter (type2 == "Coral-W")

permanova <- adonis2(OTUp_CoralW~ group ,
                     data = samples_CoralW , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_CoralW~ group ,
                                data = samples_CoralW, permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++

OTUp_MangroveW <- OTUp %>% filter (type2 == "Mangrove: Water") %>% select( -c(Row.names, type2, sample, group ))
samples_MangroveW  <- samples_dfp %>% filter (type2 == "Mangrove: Water")

permanova <- adonis2(OTUp_MangroveW ~ group ,
                     data = samples_MangroveW  , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_MangroveW~ group ,
                                data = samples_MangroveW, permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++

OTUp_MangroveS <- OTUp %>% filter (type2 == "Mangrove: Sediment") %>% select( -c(Row.names, type2, sample, group ))
samples_MangroveS <- samples_dfp %>% filter (type2 == "Mangrove: Sediment")

permanova <- adonis2(OTUp_MangroveS~ group ,
                     data = samples_MangroveS , permutations=999, method = "bray")
permanova

permanova_p <- pairwise.adonis2(OTUp_MangroveS~ group ,
                                data = samples_MangroveS, permutations=999, method = "bray")
permanova_p

############################


