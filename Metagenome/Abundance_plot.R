library(tidyverse)
library(readxl)     
library(scales)
rm(list = ls())

abund <- read_excel("Abund_species.allfilter.Metagenome.xlsx", sheet ="abundance")
samples_df <- read_excel("Abund_species.allfilter.Metagenome.xlsx", sheet ="Sample")

samples_df <- samples_df %>%  filter(!sample%in%c("qiagen_PC1", "titanF_PC3","titanL_PC2")) 
oldnames <- samples_df$F
newnames <- samples_df$sample
data.table::setnames(abund,oldnames,newnames, skip_absent=TRUE)

abunddat <- abund %>% select(Superkingdom, starts_with("qiagen"),starts_with("titanF"),starts_with("titanL")) %>%
  pivot_longer(!Superkingdom, names_to = "sample", values_to = "count")%>%
  group_by(Superkingdom, sample) %>% 
  summarise(count = sum(count))%>%
  right_join(samples_df)

source("functions.R")

col <-c( "#E6AB02","#66A61E", "#56B4E9", "#CC79A7","#596A98")

abunddat <- abunddat %>% select(Superkingdom, sample,count,  label1, type2, environment, method ) 

abund_plot (abunddat,x.var="label1",y.var="count",fill.var="Superkingdom",x=0, y=16000000, col)

pdf("Abundance_plot.pdf", width=8, height=10)
abund_plot (abunddat,x.var="label1",y.var="count",fill.var="Superkingdom",x=0, y=16000000, col)
dev.off()  
  
abunddat0 <- abunddat %>% filter(!Superkingdom == "Unclassified") 
pdf("Abundance_plot_ABVE.pdf", width=8, height=10)
abund_plot (abunddat0,x.var="label1",y.var="count",fill.var="Superkingdom",x=0, y=6000000, col)
dev.off() 

abunddat1 <- abunddat %>% filter(!Superkingdom %in% c( "Unclassified", "Eukaryota") )
col <-c( "#E6AB02","#66A61E", "#596A98")
pdf("Abundance_plot_ABV.pdf", width=8, height=10)
abund_plot (abunddat1,x.var="label1",y.var="count",fill.var="Superkingdom",x=0, y=6000000, col)
dev.off() 
#############################################################
# Relative Abundance plots
col <-c( "#E6AB02","#66A61E", "#56B4E9", "#CC79A7","#596A98")
abunddat<- abunddat %>% 
  group_by(sample)%>%
  mutate(rabund= count/sum(count))%>% 
  arrange(sample)

abund_plot (abunddat, x.var="label1",y.var="rabund",fill.var="Superkingdom",x=0, y=1, col)

pdf("Relative_abundance_kingdom.pdf", width=8, height=10)
abund_plot (abunddat, x.var="label1",y.var="rabund",fill.var="Superkingdom",x=0, y=1, col)
dev.off() 
# Relative Abundance plot excluding `Unclassified`
dat <- abunddat %>% select(Superkingdom, sample,count,  label1, type2, environment, method ) %>%
  filter(Superkingdom != "Unclassified") %>%
  group_by(sample)%>%
  mutate(relAbund= count/sum(count))%>% 
  arrange(sample)

pdf("Relative_abundance_kingdom_ABEV.pdf", width=8, height=10)
abund_plot (dat, x.var="label1",y.var="relAbund",fill.var="Superkingdom",x=0, y=1, col)
dev.off()

