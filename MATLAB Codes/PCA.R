################################################################################
# Load packages
################################################################################

library('Hmisc')
library('igraph')
library('Matrix')
library('tidyverse')
library('vegan')
library('reshape2')
library('dplyr')
library('evobiR')
library('slider')
library('runner')
library('readxl')
library('writexl')
library('xlsx')
library(ggbiplot)
library(pracma)
library(plyr)
library(gridExtra)
library(factoextra)
library(pca3d)

################################################################################
# Performing PCA analysis - in vitro + in vivo
################################################################################

# Load data
Dat <- read.csv(file = "C:/Users/fahamed2/OneDrive - University of Nebraska-Lincoln/Firnaaz Ahamed/Ongoing Research/WHONDRS Crowdsourced Manuscript/Topic6_Git/MATLAB Codes/PCAdata.csv")
grouping <- read.csv(file = "C:/Users/fahamed2/OneDrive - University of Nebraska-Lincoln/Firnaaz Ahamed/Ongoing Research/WHONDRS Crowdsourced Manuscript/Topic6_Git/MATLAB Codes/PCAgrouping.csv")

#idxRemoveSamp1 <- which(colnames(Dat)=="S19S_0079_Sed_Field_ICR_D_p2")
idxRemoveSamp <- which(grouping[,2]=="na")
#idxRemoveSamp <- unique(c(idxRemoveSamp1,idxRemoveSamp2))

Dat <- Dat[,-idxRemoveSamp]
grouping <- grouping[-idxRemoveSamp,]

# grouping
grouping1 <- grouping$grouping1
grouping2 <- grouping$grouping2
grouping3 <- grouping$grouping3

mat <- t(Dat)

# pca
idx0 <-  which(colSums(mat) %in% 0)
idxVar0 <- which(apply(mat, 2, var)==0)
idx <- unique(c(idx0,idxVar0))

if (isempty(idx)) {
} else {
  mat <- mat[,-idx]
}

mat.pca <- prcomp(mat, center=TRUE, scale. = TRUE)

# plot  
fviz_eig(mat.pca)

fviz_pca_ind(mat.pca,
             col.ind = grouping3, # color by groups
             palette = c("#00AFBB",  "#FC4E07", "#BD33A4", "#C84343"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             repel = TRUE,
             geom="point"
             )

# pca3d(mat.pca, group=grouping3)

