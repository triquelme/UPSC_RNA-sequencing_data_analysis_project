#' ---
#' title: "Heatmap"
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

#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Source multidensity function
source("~/Git/UPSCb/src/R/plot.multidensity.R")

##### vst_blind

#' load samples infos and normalized (variance stabilized) gene expression data stored in rda file
load("vst_blind.rda")

#' store the dates as factor to use as heatmap legend 
dates <- ordered_dates
dates[duplicated(dates)] <- ""
factor_dates <- factor(dates)
factor_dates
#levels(factor_dates) <- c(" ","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")

#' 1000 most variable genes showing per gene expression z-scores.
sel_var <- order(apply(vst_blind,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst_blind[sel_var,],dendrogram = "row",Colv=FALSE,labRow = NA,trace = "none",cexCol = 0.6,scale = "row",
          labCol = factor_dates)

#' ### Find the threshold for saturated expression

#' Select genes that are expressed above a treshold of counts 'exp' and in at least 'n' replicates
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

#' find the right threshold to get around 20000 genes
sapply(X = 1:10, function(i){sum(geneSelect(vst_blind,ordered_dates,i))})

#' selection of genes whose counts are over or equal to a vst value of 5 in at least 2 replicates of one time point
catalog.sel <- geneSelect(vst_blind,ordered_dates,5)
sum(catalog.sel)

#' plot cumulative gene coverage
plot(density(rowMeans(vst_blind)),col=pal[1],
     main="gene mean vst counts distribution",
     xlab="mean raw counts")

#' plot samples' variance stabilized counts distribution
cols <- sample(pal,30,replace = TRUE)
plot.multidensity(mclapply(1:ncol(vst_blind),function(k){vst_blind[catalog.sel,k]},mc.cores=16L),
                  col=cols[as.integer(factor(ordered_dates))],
                  legend.x="topright",
                  legend=levels(factor(ordered_dates)),
                  legend.col=cols,
                  legend.lwd=2,
                  main="gene vst counts distribution per sample",
                  xlab="vst")

abline(v=c(3,8))
boxplot(vst_blind)

#' Compute Saturated expression
vst.sat <- vst_blind[catalog.sel,]
vst.sat[vst.sat > 8 ] <- 8
vst.sat[vst.sat < 3 ] <- 3

#' HeatmapColourPalette to color the expression
hpal <- colorRampPalette(colors = c("blue","white","red"))(100)
# ColumnColourPalette to color the samples
cpal <- rainbow(30)

#' clusterize genes
#' with row dendrogram and without column dendrogram
heatmap.2(as.matrix(vst.sat),
          labRow=NA,trace="none",
          Colv = FALSE, dendrogram = "row",
          las=2,col=hpal,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          labCol = factor_dates,
          key=TRUE,main=" Saturated expression across chronological samples ",
          cexCol = 0.8,srtCol = 45)

#' clusterize genes and samples
#' with both row and col dendrograms
heatmap.2(as.matrix(vst.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          labCol = factor_dates,
          key=TRUE,main="Saturated expression across clusterized samples",
          cexCol = 0.8,srtCol = 45)

# NB "las=2": legend orientation vertical

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```