#' ---
#' title: "Spruce needles samples hierarchical clustering"
#' author: "Nicolas Delhomme, Thomas Riquelme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/count_data/")

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/")
#' ```
#' Libraries
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(dendextend))

#' Sources
source("~/Git/UPSCb/src/R/expressionSpecificityUtility.R")
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Plotting param
mar <- par("mar")
pal <- brewer.pal(8,"Dark2")

#' Read data
load(file = "vst_blind.rda")

#' # Prepare data

#' Replace sample id by sample dates
#' 
#' Aim: to recognize the dates and the seasons for clustering  
stopifnot(all(colnames(vst_blind) == names(ordered_dates)))
colnames(vst_blind) <- ordered_dates

#' # Clustering 
#'
#' Compute the distance of counts data between samples 
d <- dist(t(vst_blind))
# dist function computes the distances between the rows of a matrix
# we want to compute the distances between samples so we transpose the count table 
# to have samples as rows

#' Compute hierarchical clusters
hc <- hclust(d)
# hclust function takes as argument a distance matrix (dissimilarity object)

#' Plot the dendrogram
plot(hc, main = "Hierarchical clustering of samples labeled with sampling dates")


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```