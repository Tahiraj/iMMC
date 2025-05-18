library(tidyverse)
library(phyloseq)
library(gridExtra)
library(vegan)
library(pairwiseAdonis)

rm(list = ls())

immc_decontam <- readRDS("immc_decontam.rds")
source("functions.R")
################## Coral
sn = sample_names(sample_data(immc_decontam)[sample_data(immc_decontam)$type2=="Coral1" |sample_data(immc_decontam)$type2=="Coral2"|sample_data(immc_decontam)$type2=="Coral-W" , ])

rb = prune_samples(sample_names(immc_decontam) %in% sn, immc_decontam)

p1 <- betadiv_plot(rb, 1, scale="fixed") 
################## Mangrove
sn = sample_names(sample_data(immc_decontam)[sample_data(immc_decontam)$type2=="Mangrove: Water" |sample_data(immc_decontam)$type2=="Mangrove: Sediment", ])

rb = prune_samples(sample_names(immc_decontam) %in% sn, immc_decontam)

p2 <- betadiv_plot(rb, 1, scale="fixed") 

#pdf("Beta_decontam_plot.pdf", width=12, height=8)
grid.arrange(p1 ,p2, nrow=2, widths=c(1,  1, 0.6 ), layout_matrix = rbind(c(1,1,1), c(2, 2,NA)))
#dev.off()
############################################################################################# 
# Pemanova test

OTU_decontaminated <- otu_table(immc_decontam)
OTUp <- t(OTU_decontaminated) 

samples_df <- data.frame(sample_data(immc_decontam))  
samples_dfp <- samples_df %>% select(type2, sample, group )

OTUp <- merge(samples_dfp , OTUp, by = 'row.names', all = TRUE) 

#++++++++++++++++++++++++++++++ Coral1
OTUp_Coral1 <- OTUp %>% filter (type2 == "Coral1") %>% select( -c(Row.names, type2, sample, group ))
samples_Coral1 <- samples_dfp %>% filter (type2 == "Coral1")

permanova <- adonis2(OTUp_Coral1~ group ,
                     data = samples_Coral1 , permutations=999, method = "bray")
permanova
permanova_p <- pairwise.adonis2(OTUp_Coral1~ group ,
                                data = samples_Coral1 , permutations=999, method = "bray")
permanova_p
#++++++++++++++++++++++++++++++ Coral2

OTUp_Coral2 <- OTUp %>% filter (type2 == "Coral2") %>% select( -c(Row.names, type2, sample, group ))
samples_Coral2 <- samples_dfp %>% filter (type2 == "Coral2")

permanova <- adonis2(OTUp_Coral2~ group ,
                     data = samples_Coral2 , permutations=999, method = "bray")

permanova
permanova_p <- pairwise.adonis2(OTUp_Coral2~ group ,
                                data = samples_Coral2 , permutations=999, method = "bray")
permanova_p
#++++++++++++++++++++++++++++++ Coral-Water

OTUp_CoralW <- OTUp %>% filter (type2 == "Coral-W") %>% select( -c(Row.names, type2, sample, group ))
samples_CoralW <- samples_dfp %>% filter (type2 == "Coral-W")

permanova <- adonis2(OTUp_CoralW~ group ,
                     data = samples_CoralW , permutations=999, method = "bray")
permanova
permanova_p <- pairwise.adonis2(OTUp_CoralW~ group ,
                                data = samples_CoralW, permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++ Mangrove-Water

OTUp_MangroveW <- OTUp %>% filter (type2 == "Mangrove: Water") %>% select( -c(Row.names, type2, sample, group ))
samples_MangroveW  <- samples_dfp %>% filter (type2 == "Mangrove: Water")

permanova <- adonis2(OTUp_MangroveW ~ group ,
                     data = samples_MangroveW  , permutations=999, method = "bray")
permanova
permanova_p <- pairwise.adonis2(OTUp_MangroveW~ group ,
                                data = samples_MangroveW, permutations=999, method = "bray")
permanova_p

#++++++++++++++++++++++++++++++ Mangrove-Sediment

OTUp_MangroveS <- OTUp %>% filter (type2 == "Mangrove: Sediment") %>% select( -c(Row.names, type2, sample, group ))
samples_MangroveS <- samples_dfp %>% filter (type2 == "Mangrove: Sediment")

permanova <- adonis2(OTUp_MangroveS~ group ,
                     data = samples_MangroveS , permutations=999, method = "bray")
permanova
permanova_p <- pairwise.adonis2(OTUp_MangroveS~ group ,
                                data = samples_MangroveS, permutations=999, method = "bray")
permanova_p

############### Estimate_richness (Alpha Diversity)

#pdf("Diversity_plots.pdf", width=15, height=7)
plotrichness (immc_decontam, measure= "Shannon" ) 
plotrichness (immc_decontam, measure= "InvSimpson")
plotrichness (immc_decontam, measure= "Observed")
plotrichness (immc_decontam, measure= "Chao1")
#dev.off()

################ Statistical Analysis of Alpha Diversity

#Statistical Analysis of Alpha Diversity
#To estimate alpha diversity, we first evaluated the normality of the data using the Shapiro-Wilk test. 
#Significant p-values for both the inverse Simpson index and Observed ASV richness indicated a rejection of the null hypothesis, 
#confirming that the data were not normally distributed. 
#Consequently, we employed the Kruskal-Wallis test, a non-parametric alternative to ANOVA, to assess statistical significance.

diversity <- estimate_richness(immc_decontam,
                               measures = c("Observed", "Shannon"))

sample_div <- sample_data(diversity)

immc_decontam_div <- merge_phyloseq(immc_decontam, sample_div)
# Run Shapiro test
shapiro_test_Shannon <- shapiro.test(sample_data(immc_decontam_div)$Shannon)
shapiro_test_Observed <- shapiro.test(sample_data(immc_decontam_div)$Observed)

sampledataDF <- data.frame(sample_data(immc_decontam_div))

kruskal.test(Shannon ~ type2, data = sampledataDF)

#################
#Shannon
################ Mangrove-Water
sampledataDF_MW <- sampledataDF %>% filter (type2 == "Mangrove: Water") 
kruskal.test.shannon_MW <- kruskal.test(Shannon ~ group, data = sampledataDF_MW)
kruskal.test.shannon_MW
pairwise.wilcox.test(sampledataDF_MW $Shannon, sampledataDF_MW$group,
                     p.adjust.method = "fdr")

################ Mangrove-Sediment
sampledataDF_MS <- sampledataDF %>% filter (type2 == "Mangrove: Sediment") 
kruskal.test.shannon_MS <- kruskal.test(Shannon ~ group, data = sampledataDF_MS)
kruskal.test.shannon_MS
pairwise.wilcox.test(sampledataDF_MS $Shannon, sampledataDF_MS$group,
                     p.adjust.method = "fdr")

################ Coral1
sampledataDF_C1 <- sampledataDF %>% filter (type2 == "Coral1") 
kruskal.test.shannon_C1 <- kruskal.test(Shannon ~ group, data = sampledataDF_C1)
kruskal.test.shannon_C1 
pairwise.wilcox.test(sampledataDF_C1 $Shannon, sampledataDF_C1$group,
                     p.adjust.method = "fdr")

################ Coral2
sampledataDF_C2 <- sampledataDF %>% filter (type2 == "Coral2") 
kruskal.test.shannon_C2 <- kruskal.test(Shannon ~ group, data = sampledataDF_C2)
kruskal.test.shannon_C2
pairwise.wilcox.test(sampledataDF_C2 $Shannon, sampledataDF_C2$group,
                     p.adjust.method = "fdr")

################ Coral-Water
sampledataDF_CW <- sampledataDF %>% filter (type2 == "Coral-W") 
kruskal.test.shannon_CW <- kruskal.test(Shannon ~ group, data = sampledataDF_CW)
kruskal.test.shannon_CW
pairwise.wilcox.test(sampledataDF_CW $Shannon, sampledataDF_CW$group,
                     p.adjust.method = "fdr")
#################
#+Observed
################# Mangrove-Water
sampledataDF_MW <- sampledataDF %>% filter (type2 == "Mangrove: Water") 
kruskal.test.shannon_MW <- kruskal.test(Observed ~ group, data = sampledataDF_MW)
kruskal.test.shannon_MW
pairwise.wilcox.test(sampledataDF_MW $Observed, sampledataDF_MW$group,
                     p.adjust.method = "fdr")

################# Mangrove-Sediment
sampledataDF_MS <- sampledataDF %>% filter (type2 == "Mangrove: Sediment") 
kruskal.test.shannon_MS <- kruskal.test(Observed ~ group, data = sampledataDF_MS)
kruskal.test.shannon_MS
pairwise.wilcox.test(sampledataDF_MS $Observed, sampledataDF_MS$group,
                     p.adjust.method = "fdr")

################# Coral1
sampledataDF_C1 <- sampledataDF %>% filter (type2 == "Coral1") 
kruskal.test.shannon_C1 <- kruskal.test(Observed ~ group, data = sampledataDF_C1)
kruskal.test.shannon_C1
pairwise.wilcox.test(sampledataDF_C1 $Observed, sampledataDF_C1$group,
                     p.adjust.method = "fdr")

################# Coral2
sampledataDF_C2 <- sampledataDF %>% filter (type2 == "Coral2") 
kruskal.test.shannon_C2 <- kruskal.test(Observed ~ group, data = sampledataDF_C2)
kruskal.test.shannon_C2 
pairwise.wilcox.test(sampledataDF_C2 $Observed, sampledataDF_C2$group,
                     p.adjust.method = "fdr")

################# Coral-Water
sampledataDF_CW <- sampledataDF %>% filter (type2 == "Coral-W") 
kruskal.test.shannon_CW <- kruskal.test(Observed ~ group, data = sampledataDF_CW)
kruskal.test.shannon_CW
pairwise.wilcox.test(sampledataDF_CW $Observed, sampledataDF_CW$group,
                     p.adjust.method = "fdr")
#==================================================================
# Here the shapiro.test test is significant, so anvoa is not needed on shannon index
#=================================================
aov.shannon <- aov(Shannon ~ type2, data = sampledataDF)
summary(aov.shannon)
TukeyHSD(aov.shannon)$type2
#+++++++++++++++++++++++++++++++++++++++++
sampledataDF_MW <- sampledataDF %>% filter (type2 == "Mangrove: Water") 
aov.shannon_MW <- aov(Shannon ~ group, data = sampledataDF_MW)
summary(aov.shannon_MW)
TukeyHSD(aov.shannon_MW)$group

