library(tidyverse)
library(phyloseq)
library(readxl)       #to import the data from Excel file
library(scales)
rm(list = ls())

abund <- read_excel("abundance_table_species_ncbi_q15.xlsx", sheet ="abundance")
samples_df <- read_excel("abundance_table_species_ncbi_q15.xlsx", sheet ="Sample")

abund_Archaea <- abund %>% select(Superkingdom,starts_with("qiagen"),starts_with("titanF"),starts_with("titanL")) %>%
  pivot_longer(!Superkingdom, names_to = "sample", values_to = "count")%>%
  group_by(Superkingdom, sample) %>% 
  summarise(Archaea_count = sum(count))%>%
  right_join(samples_df) %>% filter (Superkingdom=="Archaea" & Archaea_count>0) %>%
  select(Superkingdom, sample, Archaea_count, label1)

abunddat <- abund %>% select(Superkingdom, starts_with("qiagen"),starts_with("titanF"),starts_with("titanL")) %>%
  pivot_longer(!Superkingdom, names_to = "sample", values_to = "count")%>%
  group_by(Superkingdom, sample) %>% 
  summarise(count = sum(count))%>%
  right_join(samples_df) %>% 
  left_join(abund_Archaea %>% select(sample, Archaea_count)) %>% 
  filter (count > 0)

Abundance_plot <- function(dat, labels, limit, cols) {
  p= ggplot(dat, aes(x=label1, y= count, fill=Superkingdom)) +
    geom_point(aes( y=1200000, x=label1, size=Archaea_count), col="#E6AB02")+
    geom_bar(stat="identity", position="stack", width=0.5)+
    # geom_bar(stat="identity", position=position_dodge())+
    scale_fill_manual("",values = cols)+
    theme_minimal() + xlab("")+ ylab("Abundance") +
    scale_y_continuous(" ", limits = limit, n.breaks =6,labels = scales::comma)+
    guides(size = guide_legend(title ="Archaea"))+
    facet_grid(~factor(method, levels=labels) ,   scales = "free_x")
    #      labels = c("Coral" /n "Qiagen (Lab)", "Mangrove" /n "Qiagen (Lab))")),
  p
}

cols <- c("#E6AB02","#66A61E", "#56B4E9", "#CC79A7")
dat <-  abunddat %>% filter (environment=="Lab_Q") 
p1 <- Abundance_plot (dat = dat, labels = c("Coral: Qiagen (Lab)", "Mangrove: Qiagen (Lab)") , 
                      limit = c(0, 1250000), cols)
p1 <- p1 + theme(legend.position = "none", axis.text.x=element_blank())

dat <-  abunddat %>% filter (environment=="Field") 
p2 <- Abundance_plot (dat = dat, labels = c("xTitan (Boat)" , "xTitan (In situ)"), 
                      limit = c(0, 1250000), cols)
p2 <- p2 + theme(legend.position = "none", axis.text.x=element_blank())

dat <-  abunddat %>% filter (environment=="Lab") 
p3 <- Abundance_plot (dat = dat, labels = c("xTitan (Lab)" , "X-Titan (Lab)"), 
                      limit = c(0, 1250000), cols) 
p3 <- p3 + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("Abundance_plot.pdf", width=8, height=10)
gridExtra::grid.arrange(p1,p2,p3, heights=c(1,1,1.5))
dev.off() 
#############################################################
# Log-Abundance filter the no Eukaryota

abunddat0 <- abunddat %>% filter(Superkingdom !="Eukaryota") 
 
Abundance_plot_log <- function(dat, labels, limit, cols) {
  
 p <- ggplot(dat , aes(x=label1, y= count, fill=factor(Superkingdom, levels=c("Bacteria","Archaea", "Unclassified")))) +
  geom_bar(stat="identity", position="stack", width=0.5)+
  scale_fill_manual("",values = cols)+
  theme_minimal() + xlab("")+ ylab("Abundance") +
  scale_y_log10(" ",limits = limit, breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(~factor(method, levels=labels) ,   scales = "free_x")

}

cols <- c("#66A61E", "#E6AB02", "#CC79A7")

dat <-  abunddat0 %>% filter (environment=="Lab_Q") 
p1 <- Abundance_plot_log (dat, labels = c("Coral: Qiagen (Lab)", "Mangrove: Qiagen (Lab)") , 
                          limit = c(1,1e11), cols )
p1 <- p1 + theme(legend.position = "none", axis.text.x=element_blank())

dat <-  abunddat0 %>% filter (environment=="Field") 
p2 <- Abundance_plot_log (dat, labels = c("xTitan (Boat)" , "xTitan (In situ)"), 
                          limit = c(1,1e11), cols  )
p2 <- p2 + theme(legend.position = "none", axis.text.x=element_blank())

dat <-  abunddat0 %>% filter (environment=="Lab") 
p3 <- Abundance_plot_log (dat, labels = c("xTitan (Lab)" , "X-Titan (Lab)"), 
                          limit = c(1,1e11), cols) 
p3 <- p3 + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gridExtra::grid.arrange(p1,p2,p3, heights=c(1,1,1.5))

#####################################################

abund_per <- abund %>% mutate(
  Qiagen = select(., starts_with("qiagen")) %>% rowSums(),
  TitanF = select(., starts_with("titanF")) %>% rowSums(),
  TitanL = select(., starts_with("titanL")) %>% rowSums()) %>%
  select(Superkingdom, Qiagen, TitanF, TitanL )%>%
  pivot_longer(!Superkingdom, names_to = "method", values_to = "count")%>%
  group_by(Superkingdom, method) %>% 
  summarise(count = sum(count)) %>% 
  group_by(method) %>% 
  mutate(perc = 100*count/sum(count))%>%
  ungroup() 

ggplot(abund_per, aes(x=method, y= perc, fill=Superkingdom)) +
  geom_bar(stat="identity", position="stack", width=0.5)+
  scale_fill_manual("",values = c("#F0E442","#66A61E", "#56B4E9", "#CC79A7"))+
  theme_minimal() + xlab("")+ ylab("Abundance (%)")

#===========================================================================


