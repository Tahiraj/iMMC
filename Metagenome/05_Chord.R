
library(circlize)
library(ComplexHeatmap)

rm(list = ls())
source("functions.R")
dat <- read_csv("dat_chord_mangrove.csv")

dat %>% pivot_longer(cols= Aliifodinibius:Fodinibius, names_to = "Genus", values_to = "count")

#pdf("chord.pdf", width=8, height=8)
###################### Mangrove
orders <- c("Qiagen", "titanF", "titanL",
            "Fodinibius", "Rhodohalobacter", "Gracilimonas","Aliifodinibius", "Balneola", "Unclassified Balneolaceae")

col = c("darkcyan", "darkgoldenrod2", "#FC4E07")
####################### Mangrove: Sediment

dat1 <- dat %>% filter(Biome == "Sediment")
mat1 <- dat1 %>% select(-c(Biome, Method ))%>%
  select_if(~ sum(.x) > 1) %>% as.matrix()

rownames(mat1) <- dat1$Method

col_mangrove = c("#FFE528", "#8E0152", "#EBA3CC","#7570B3",  "#D3ECB2", "green4")

chord(mat1, 
      grid.col=c(col, col_mangrove),
      orders, 
      col, label = FALSE)

LM = Legend(title = "Mangrove-Sediment", 
            labels = c("Fodinibius","Rhodohalobacter", "Gracilimonas","Aliifodinibius", "Balneola", "Unclassified Balneolaceae"),
            legend_gp = gpar(fill=col_mangrove,
                             col = "white", cex = 0.7))
draw(LM, x = unit(2, "cm"), y = unit(2, "cm"), just = c("left", "bottom"))

####################### Mangrove: Seawater

dat2 <- dat %>% filter(Biome == "Seawater")
mat2 <- dat2 %>% select(-c(Biome, Method ))%>%
  select_if(~ sum(.x) > 1) %>% as.matrix()

rownames(mat2) <- dat2$Method
chord(mat2, 
      grid.col=c(col, col_mangrove),
      orders, 
      col, label = FALSE)
LSW = Legend(title = "Mangrove-Water", 
             labels = c("Fodinibius","Rhodohalobacter", "Gracilimonas","Aliifodinibius", "Balneola", "Unclassified Balneolaceae"),
            legend_gp = gpar(fill=col_mangrove,
                             col = "white", cex = 0.7))
draw(LSW, x = unit(2, "cm"), y = unit(2, "cm"), just = c("left", "bottom"))

########################## Coral

dat_Coral <- read_csv("dat_chord_coral.csv")

###################### Coral1
dat_Coral1 <- dat_Coral %>% filter(Biome == "Coral1")
matcoral1 <- dat_Coral1 %>% select(-c(Biome, Method ))%>%
  select_if(~ sum(.x) > 1) %>% as.matrix()

rownames(matcoral1) <- dat_Coral1 $Method

col_coral1 <- c("#8E0152", "#EBA3CC", "#7570B3", "#D3ECB2", "green4","#FFE528")

chord(matcoral1, 
      grid.col=c( col, col_coral1), 
      orders=c(rownames(matcoral1), colnames(matcoral1)), 
      col, label = FALSE)

L1 = Legend(title = "Coral1", labels = colnames(matcoral1),
            legend_gp = gpar(fill=col_coral1,
                             col = "white", cex = 0.7))
draw(L1, x = unit(2, "cm"), y = unit(2, "cm"), just = c("left", "bottom"))

###################### Coral2. 
dat_Coral2 <- dat_Coral %>% filter(Biome == "Coral2")
matcoral2 <- dat_Coral2 %>% select(-c(Biome, Method ))%>%
  select_if(~ sum(.x) > 1) %>% as.matrix()
rownames(matcoral2) <- dat_Coral2 $Method

col_coral2 = c( "#8E0152", "#EBA3CC", "#7570B3", "#D3ECB2", "green4", "#FFE528", "magenta", "lightpink4")

chord(matcoral2, 
      grid.col=c( col, col_coral2), 
      orders= c(rownames(matcoral2), colnames(matcoral2)), 
      col, label = FALSE)
L2 = Legend(title = "Coral2", labels = colnames(matcoral2),
            legend_gp = gpar(fill=col_coral2,
                             col = "white", cex = 0.7))
draw(L2, x = unit(2, "cm"), y = unit(2, "cm"), just = c("left", "bottom"))

###################### Coral Water

dat_Coralw <- dat_Coral %>% filter(Biome == "Coral-W")
matcoralw <- dat_Coralw %>% select(-c(Biome, Method ))%>%
  select_if(~ sum(.x) > 1) %>% as.matrix()
rownames(matcoralw) <- dat_Coralw $Method

col_coralw = c( "#8E0152", "#EBA3CC", "#7570B3", "#D3ECB2", "green4", "#FFE528", "magenta")

chord(matcoralw, 
      grid.col=c( col, col_coralw), 
      orders=c(rownames(matcoralw), colnames(matcoralw)), 
      col, label = FALSE)

L3 = Legend(title = "Coral-Water", labels = colnames(matcoral1),
            legend_gp = gpar(fill=col_coralw,
                             col = "white", cex = 0.7))
draw(L3, x = unit(2, "cm"), y = unit(2, "cm"), just = c("left", "bottom"))

####################

pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect()  # border
draw(L1, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
draw(L2, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
draw(L3, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
popViewport()
