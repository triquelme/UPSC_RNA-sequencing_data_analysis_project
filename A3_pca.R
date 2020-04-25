#' ---
#' title: "Principal Component Analysis"
#' author: "Nicolas Delhomme & Thomas Riquelme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/count_data/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```
#' Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
library(tidyverse)
library(plotly)

###' PCA for vst blind counts
#' load samples infos and normalized (variance stabilized) gene expression data stored in rda file
load("vst_blind.rda")

#' Compute principal component analysis
pc <- prcomp(t(vst_blind)) 
# 't' transposes a matrix which means switch rows with columns
# before in vst.kt: rows=genes and columns=samples and data inside  = counts                    
# here PCA is calculated row after row, and we want to calculate it for each sample (so 84 PCs)

percent <- round(summary(pc)$importance[2,]*100)
#take the proportion of variance *100 and round for every PC
percent

#' plot according to PC1 and PC2  
p <- data.frame(sample = rownames(pc$x), pc$x[, 1:2], date = ordered_dates) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = date)) +
  theme(legend.position = "right") +
  geom_point()

ggplotly(p)

ggsave("pca_plot2_PC1_PC2.png")

#' plot according to PC2 and PC3

p2 <- data.frame(sample = rownames(pc$x), pc$x[, 2:3], 
                 date = ordered_dates) %>% 
  ggplot(aes(x = PC2, y = PC3, colour = date)) +
  theme(legend.position = "right") +
  geom_point()

ggplotly(p2)

ggsave("pca_plot2_PC2_PC3.png")


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```