library(tidyverse)
library(phyloseq)
library(readxl)       
library(decontam)
#decontam package provides simple statistical methods to identify and visualize contaminating DNA features
packageVersion("phyloseq")
rm(list = ls())

abund <- read_excel("Abund_species.allfilter.Metagenome.xlsx", sheet ="abundance")
samples_dat <- read_excel("Abund_species.allfilter.Metagenome.xlsx", sheet ="Sample")%>% 
  filter(!sample%in%c("qiagen_PC1", "titanF_PC3","titanL_PC2")) 
oldnames <- samples_dat$F
newnames <- samples_dat$sample
data.table::setnames(abund,oldnames,newnames, skip_absent=TRUE)

abund_dat <- abund%>% 
  filter(Superkingdom %in% c("Archaea", "Bacteria", "Viruses") )%>%
  # filter(Superkingdom == "Bacteria") %>% 
  select(Superkingdom, Species, starts_with("qiagen"),starts_with("titanF"),starts_with("titanL"))%>%
  mutate(across(where(is.numeric), function(x) ifelse(x <= 10, 0, x)))

otu_dat <- abund_dat %>% 
  rowwise() %>% 
  filter(sum(c_across("qiagen_C13":"titanL_M36")) != 0) 

tax_dat <- read_excel("Abund_species.allfilter.Metagenome.xlsx", sheet ="taxa")
tax_dat <- tax_dat %>% filter (Species %in% otu_dat$Species)

samples_df <- samples_dat %>%
  tibble::column_to_rownames("sample") 

colnames(samples_df)[1] <- "sample"

OTU = otu_table(as.matrix(otu_dat[,-c(1:2)]), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_dat))
samples = sample_data(samples_df)

immc_ps <- phyloseq(OTU, TAX, samples)

sample_data(immc_ps)$is.neg <- sample_data(immc_ps)$type == "Neg Control"

contamdf.prev05 <- isContaminant(immc_ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
which(contamdf.prev05$contaminant)

contamdf.prev01 <- isContaminant(immc_ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
which(contamdf.prev01$contaminant)

#=================================================================
contaminated_rows <- which(contamdf.prev05$contaminant)
OTU_contaminated <- OTU[contaminated_rows,]

TAX_contaminated = TAX[contaminated_rows,]

contaminated_dat <- merge(TAX_contaminated, OTU_contaminated, by = 'row.names', all = TRUE) 

OTU_decontaminated <- OTU[-contaminated_rows,]

TAX_decontaminated = TAX[-contaminated_rows,]

#=================================================================

immc_decontam <- phyloseq(OTU_decontaminated, TAX_decontaminated, samples)

#saveRDS(immc_decontam, "immc_decontam.rds")
