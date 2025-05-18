library(tidyverse)
library(phyloseq)
library(readxl)       
library(decontam)
packageVersion("phyloseq")
#decontam package provides simple statistical methods to identify and visualize contaminating DNA features
rm(list = ls())

abund <- read_excel("abundance_table_species_ncbi_q15.xlsx", sheet ="abundance")

abund_dat <- abund%>% filter(Superkingdom == "Bacteria") %>% 
  select(Superkingdom, Species, starts_with("qiagen"),starts_with("titanF"),starts_with("titanL"))%>%
  mutate(across(where(is.numeric), function(x) ifelse(x <= 10, 0, x)))

otu_dat <- abund_dat %>% 
  rowwise() %>% 
  filter(sum(c_across(qiagen_C13:titanL_PC2)) != 0) %>%
  select(!c(qiagen_PC1, titanF_PC3, titanL_PC2))

tax_dat <- read_excel("abundance_table_species_ncbi_q15.xlsx", sheet ="taxa")
tax_dat <- tax_dat %>% filter (Species %in% otu_dat$Species)

samples_df <- read_excel("abundance_table_species_ncbi_q15.xlsx", sheet ="Sample")
samples_df <- samples_df %>% 
  filter(!sample0%in%c("qiagen_PC1", "titanF_PC3"," titanL_PC2")) %>%
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
#[1]    4   82   92  119  832  920 1021 1097 1103 1638 1962 1966 1979 2058 2068 2685 3075 4620

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
