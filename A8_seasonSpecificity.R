#' ---
#' title: "Spruce needles season specificity analysis"
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
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/")

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
load(file = "vst_aware_filtered.rda")

#' # Prepare data

#' Replace sample id by sample dates
#' 
#' Aim: to recognize the dates and the seasons for clustering  
stopifnot(all(colnames(vst_aware_filtered) == names(ordered_dates)))
colnames(vst_aware_filtered) <- ordered_dates

#' # Clustering 
#' Based on the distance of counts data between samples
hc <- hclust(dist(t(vst_aware_filtered)))
plot(hc, main = "Hierarchical clustering of samples labeled with sampling dates")

#' # Identification of season periods between clustered samples
#' Cut the tree in different clusters
ct <- cutree(hc, h = 70)
#' Reassignment of season number identifiers for each samples.
#' 
#' Because there are more samples of May and June (half of the samples), the first ramifications and clusters in the trees are between samples of these months.
#'
#' To avoid that, we merge the ramifications for these months and focus on the clusters between seasons.
ct[ct ==1 | ct == 2 | ct ==3 | ct ==4 | ct ==5 | ct ==11 | ct == 12] <- 1 # 05-06 = early-summer
ct[ct == 6] <- 2 # 07-08 = late-summer
ct[ct == 7] <- 3 # 09 = autumn
ct[ct == 8 | ct == 9] <- 4 # 10-02 = winter
ct[ct ==10] <- 5 # 03-04 == spring

#' Rename clusters with seasons'names
tp <- factor(ct,labels=c("early-summer","late-summer","autumn","winter","spring"))
tp_char <- as.character(tp)
tp_char_ordered <- tp_char[hc$order]

#' Change labels in the dendrogram
labels(hc) <- tp_char_ordered
plot(hc, main="Hierarchical clustering of samples labeled with seasons")

#' # Season expression specificity for all genes

#' ## Compute season expression specificity
es <- expressionSpecificity(vst_aware_filtered, tp_char, mode = "local", output="complete")
#' This expressionSpecificity function calculates the season expression specificity score. 
#' The score ranges from 0 to 1 i.e. from ubiquitous to season specific
#'
#' The expressionSpecificity function takes as argument an expression matrix, a factor containing the corresponding season for every sample,and two options:
#' 
#' "local" mode which computes season specificity, and  "global mode" which computes sample specificity.
#' 
#' "complete mode" returns the specificity scores for every genes, plus the average expression of a gene in each season, the maximal average expression of a gene in a season and the number of seasons taken into account for calculation where the gene expresion is higher than 0.
#'
#' In short, this function will split the samples by seasons and compute the average expression of every gene (aij) in every seaon.
#' Then, it will compute the maximal average expression in a season for every gene (max(ai in 1:n)).
#' Finally, it will calculate the score by dividing average expression in every season by the maximal average expression in a season, summing all the seasons, and then dividing by the number of season.
#'
#' ## Plot gene specificity score distribution of all genes
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
# plot margin: mar=c(bottom,left,top,right)
# axis margin: mgp=c(axis title, axis labels, axis line)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity distribution", 
        ylab = "Number of genes", xlab = "Gene season specificity score", cex.names = 0.9)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity logarithmic distribution",log="y", ylab = "Number of genes", xlab = "Season specificity score", cex.names = 0.9)
# "cut" function splits gene season specificity score values into ten levels or bins.
# Then table represents in a bar plot the number of genes that have this level of season specificity.

#' We observed that most of the genes have low or average score, which means that are very ubiquitous.
#' Only a very small fraction of them are season specific (score between 0.9 and 1), around 200 genes.
#'
#'  ## Boxplot gene expression for each gene season specificity score category
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
boxplot(split(es[,2:6],cut(es[,"score"],10)),las=2, main= "Gene average expression relatively to gene season specificity score", 
        ylab = "Gene average expression", xlab = "Gene season specificity score",cex.axis=0.9) 
# "es[,2:6]" is the average expression in individual seasons: autumn, early-summer, late-summer, spring, winter.
# Then it splits genes'average expression in individual seasons relatively to gene's specificity score level organised in 10 bins.
# Finally, it plots genes'average expression for every genes'season specificity score level.

#' Interestingly the genes with less season specificity (ie the most ubiquitous) are the most expressed.
#'
#' ## Check how many genes are never expressed in a season
pander(round(colSums(is.na(es[,2:6])) / nrow(es) * 100,2))
# "colSums(is.na(es[,2:6])" counts all the non expressed genes in a season (their average expression value is listed as na).
# The number of non expressed genes is then divided by the number of genes and multiplied by 100 to calculate the percentage.

#' Less than 5% of the genes are never expressed in any season.
#' Interestingly all of the genes are expressed during early summer. 
#' Either this could imply that all genes are necessary at that time of the year in a very active metabolism,
#' or either it could be a sampling biais as half of the samples were sampled in that time of the year.
#'
#' ## Check for genes more expressed during a season
#' For every season, select genes whose average expression in a season is equal to the maximal average expression across all seasons.
#' This list of genes obtained corresponds to the genes with an expression peak in the season tested.
#' It also plots gene season specificity distribution per season.

dev.null <- sapply(levels(tp),function(x){
  sel <- es[,x] == es[,"maxn"]
  message(sprintf("There are %s genes with an expression peak in %s",
                  sum(sel,na.rm=TRUE),x))
  par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
  barplot(table(cut(es[sel,"score"],10)),las=2,
          main=paste0(x," gene season specificity distribution"),
          ylab= "Number of genes", xlab= "Gene season specificity score", cex.names = 0.9)
})

#' ## Get lists of genes more expressed during a season
season_list <- sapply(levels(tp), function(x){
  sel <- es[,x] == es[,"maxn"]
  sel <- sel[!is.na(sel)]
  gene_list <- names(sel[sel==TRUE])
  return(gene_list)
})
par(mar=mar, mgp = c(3, 1, 0))
barplot(lengths(season_list), col = c("yellow","orange","red","blue","green"), ylab= "number of specific genes", main = "Number of specific genes per season")

#' The lists of season specific genes are too large.
#' There is a need for a threshold to focus only of the most season specific genes with the highest scores.
#'
#' Save season specific lists of genes
dir.create(file.path("analysis","season"),showWarnings = FALSE)
dev.null <- sapply(levels(tp), function(x){
  write(season_list[[x]],file = paste0("analysis/season/",x,"_peak_genes"))
})

#' ## Get lists of genes more expressed during a season with high specificity

#' Get the names of genes that are more expressed in a season than in the other, 
#' and whose score is above a threshold specificity score of 0.6.
season_specific_genes <- sapply(levels(tp), function(x){
  x.genes <- names(which(es[es[,x] == es[,"maxn"],"score"]>0.6))
  message(sprintf("There are %s genes that are %s specific",length(x.genes),x))
  return(x.genes)
})

#' Save season specific lists of genes
dev.null <- sapply(levels(tp), function(x){
  write(season_specific_genes[[x]],file = paste0("analysis/season/",x,"_specific_genes"))
})

#' # Redo season specificity analysis for photosynthetic goi

#' Read Gene Of Interest csv file
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_list.csv", header = FALSE)

#' Select counts only for those genes
sel <- vst_aware_filtered[rownames(vst_aware_filtered) %in% goi$V1,]

#' ## Compute season expression specificity
es <- expressionSpecificity(sel, tp_char, mode = "local", output="complete")

#' ## Plot gene specificity score distribution of all genes
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
# plot margin: mar=c(bottom,left,top,right)
# axis margin: mgp=c(axis title, axis labels, axis line)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity distribution", 
        ylab = "Number of genes", xlab = "Gene season specificity score", cex.names = 0.9)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity logarithmic distribution",log="y", ylab = "Number of genes", xlab = "Season specificity score", cex.names = 0.9)
# "cut" function splits gene season specificity score values into ten levels or bins.
# Then table represents in a bar plot the number of genes that have this level of season specificity.

#' We observed that most of the genes have low or average score, which means that are very ubiquitous.
#' Only a very small fraction of them are season specific (score between 0.9 and 1), around 200 genes.
#'
#'  ## Boxplot gene expression for each gene season specificity score category
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
boxplot(split(es[,2:6],cut(es[,"score"],10)),las=2, main= "Gene average expression relatively to gene season specificity score", 
        ylab = "Gene average expression", xlab = "Gene season specificity score",cex.axis=0.9) 
# "es[,2:6]" is the average expression in individual seasons: autumn, early-summer, late-summer, spring, winter.
# Then it splits genes'average expression in individual seasons relatively to gene's specificity score level organised in 10 bins.
# Finally, it plots genes'average expression for every genes'season specificity score level.

#' Interestingly the genes with less season specificity (ie the most ubiquitous) are the most expressed.
#'
#' ## Check how many genes are never expressed in a season
pander(round(colSums(is.na(es[,2:6])) / nrow(es) * 100,2))
# "colSums(is.na(es[,2:6])" counts all the non expressed genes in a season (their average expression value is listed as na).
# The number of non expressed genes is then divided by the number of genes and multiplied by 100 to calculate the percentage.

#' Less than 5% of the genes are never expressed in any season.
#' Interestingly all of the genes are expressed during early summer. 
#' Either this could imply that all genes are necessary at that time of the year in a very active metabolism,
#' or either it could be a sampling biais as half of the samples were sampled in that time of the year.
#'
#' ## Check for genes more expressed during a season
#' For every season, select genes whose average expression in a season is equal to the maximal average expression across all seasons.
#' This list of genes obtained corresponds to the genes with an expression peak in the season tested.
#' It also plots gene season specificity distribution per season.

dev.null <- sapply(levels(tp),function(x){
  sel <- es[,x] == es[,"maxn"]
  message(sprintf("There are %s genes with an expression peak in %s",
                  sum(sel,na.rm=TRUE),x))
})

#' There are 7 genes with an expression peak in early-summer
#' There are 34 genes with an expression peak in late-summer
#' There are 9 genes with an expression peak in autumn
#' There are 12 genes with an expression peak in winter
#' There are 18 genes with an expression peak in spring
#' 
#' It could be interesting to look at the photosynthetic genes that are more expressed in Autumn and Winter. 
#'
#' ## Get lists of genes more expressed during a season
season_list <- sapply(levels(tp), function(x){
  sel <- es[,x] == es[,"maxn"]
  sel <- sel[!is.na(sel)]
  gene_list <- names(sel[sel==TRUE])
  return(gene_list)
})
#' The lists of season specific genes are too large.
#' There is a need for a threshold to focus only of the most season specific genes with the highest scores.
#'
#' Save season specific lists of genes
dir.create(file.path("analysis/season/", "photosynthetic_goi"))
dev.null <- sapply(levels(tp), function(x){
  write(season_list[[x]],file = paste0("analysis/season/photosynthetic_goi/",x,"_peak_genes"))
})

#' ## Get lists of genes more expressed during a season with high specificity
season_specific_genes <- sapply(levels(tp), function(x){
  x.genes <- names(which(es[es[,x] == es[,"maxn"],"score"]>0.6))
  message(sprintf("There are %s genes that are %s specific",length(x.genes),x))
  return(x.genes)
})

#' There are 0 genes that are early-summer specific
#' There are 0 genes that are late-summer specific
#' There are 0 genes that are autumn specific
#' There are 0 genes that are winter specific
#' There are 0 genes that are spring specific
#' No genes specific enough for a season to be above the threshold of specificity score
#' 
#' # Redo season specificity analysis for cytoscape cluster 2,2 enriched in photosynthetic goi

#' Read Gene Of Interest csv file
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/cluster_goi_2,2_gene_list.csv", header = FALSE)

#' Select counts only for those genes
sel <- vst_aware_filtered[rownames(vst_aware_filtered) %in% goi$V1,]

#' ## Compute season expression specificity
es <- expressionSpecificity(sel, tp_char, mode = "local", output="complete")

#' ## Plot gene specificity score distribution of all genes
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
# plot margin: mar=c(bottom,left,top,right)
# axis margin: mgp=c(axis title, axis labels, axis line)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity distribution", 
        ylab = "Number of genes", xlab = "Gene season specificity score", cex.names = 0.9)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity logarithmic distribution",log="y", ylab = "Number of genes", xlab = "Season specificity score", cex.names = 0.9)
# "cut" function splits gene season specificity score values into ten levels or bins.
# Then table represents in a bar plot the number of genes that have this level of season specificity.

#' We observed that most of the genes have low or average score, which means that are very ubiquitous.
#' Only a very small fraction of them are season specific (score between 0.9 and 1), around 20 genes.
#'
#'  ## Boxplot gene expression for each gene season specificity score category
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
boxplot(split(es[,2:6],cut(es[,"score"],10)),las=2, main= "Gene average expression relatively to gene season specificity score", 
        ylab = "Gene average expression", xlab = "Gene season specificity score",cex.axis=0.9) 
# "es[,2:6]" is the average expression in individual seasons: autumn, early-summer, late-summer, spring, winter.
# Then it splits genes'average expression in individual seasons relatively to gene's specificity score level organised in 10 bins.
# Finally, it plots genes'average expression for every genes'season specificity score level.

#' Interestingly the genes with less season specificity (ie the most ubiquitous) are the most expressed.
#'
#' ## Check how many genes are never expressed in a season
pander(round(colSums(is.na(es[,2:6])) / nrow(es) * 100,2))
# "colSums(is.na(es[,2:6])" counts all the non expressed genes in a season (their average expression value is listed as na).
# The number of non expressed genes is then divided by the number of genes and multiplied by 100 to calculate the percentage.

#' 30% of the genes are never expressed in autumn and spring. 
#' Those seasons are also the shorter one (composed of 1 or 2 months).
#' Interestingly all of the genes are expressed during summer and almost all (96%) are expressed during winter.
#' Those seasons are also the longer ones (composed of 4 or 5 months).
#' This could be a biais introduced by the length of the seasons choosed subjectively 
#' by looking at how the samples at different times of the year cluster together in term of genes expression.
#' 
#' ## Check for genes more expressed during a season
#' For every season, select genes whose average expression in a season is equal to the maximal average expression across all seasons.
#' This list of genes obtained corresponds to the genes with an expression peak in the season tested.
#' It also plots gene season specificity distribution per season.

dev.null <- sapply(levels(tp),function(x){
  sel <- es[,x] == es[,"maxn"]
  message(sprintf("There are %s genes with an expression peak in %s",
                  sum(sel,na.rm=TRUE),x))
})

#' There are 152 genes with an expression peak in early-summer
#' There are 75 genes with an expression peak in late-summer
#' There are 13 genes with an expression peak in autumn
#' There are 59 genes with an expression peak in winter
#' There are 14 genes with an expression peak in spring
#' 
#' It could be interesting to look at the photosynthetic genes that are more expressed in Autumn and Winter.
#'
#' ## Get lists of genes more expressed during a season
season_list <- sapply(levels(tp), function(x){
  sel <- es[,x] == es[,"maxn"]
  sel <- sel[!is.na(sel)]
  gene_list <- names(sel[sel==TRUE])
  return(gene_list)
})

#' Save season specific lists of genes
dir.create(file.path("analysis/season/", "cytoscape_cluster_photosynthetic_genes_enriched_1"))
dev.null <- sapply(levels(tp), function(x){
  write(season_list[[x]],file = paste0("analysis/season/cytoscape_cluster_photosynthetic_genes_enriched_1/",x,"_peak_genes"))
})

#' The lists of season specific genes are too large.
#' There is a need for a threshold to focus only of the most season specific genes with the highest scores.
#'
#' ## Get lists of genes more expressed during a season with high specificity

#' Get the names of genes that are more expressed in a season than in the other, 
#' and whose score is above a threshold specificity score of 0.6.
season_specific_genes <- sapply(levels(tp), function(x){
  x.genes <- names(which(es[es[,x] == es[,"maxn"],"score"]>0.6))
  message(sprintf("There are %s genes that are %s specific",length(x.genes),x))
  return(x.genes)
})

#' There are 35 genes that are early-summer specific
#' There are 2 genes that are late-summer specific
#' There are 0 genes that are autumn specific
#' There are 5 genes that are winter specific
#' There are 0 genes that are spring specific
#' 
#' It could be interesting to look at the five photosynthetic genes highly specific for winter

#' Save season specific lists of genes
dev.null <- sapply(levels(tp), function(x){
  write(season_specific_genes[[x]],file = paste0("analysis/season/cytoscape_cluster_photosynthetic_genes_enriched_1/",x,"_specific_genes"))
})

#' # Redo season specificity analysis for cytoscape cluster 2,4 enriched in photosynthetic goi

#' Read Gene Of Interest csv file
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/cluster_goi_2,4_gene_list.csv", header = FALSE)

#' Select counts only for those genes
sel <- vst_aware_filtered[rownames(vst_aware_filtered) %in% goi$V1,]

#' ## Compute season expression specificity
es <- expressionSpecificity(sel, tp_char, mode = "local", output="complete")

#' ## Plot gene specificity score distribution of all genes
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
# plot margin: mar=c(bottom,left,top,right)
# axis margin: mgp=c(axis title, axis labels, axis line)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity distribution", 
        ylab = "Number of genes", xlab = "Gene season specificity score", cex.names = 0.9)
barplot(table(cut(es[,"score"],10)),las=2,main="Gene season specificity logarithmic distribution",log="y", ylab = "Number of genes", xlab = "Season specificity score", cex.names = 0.9)
# "cut" function splits gene season specificity score values into ten levels or bins.
# Then table represents in a bar plot the number of genes that have this level of season specificity.

#' We observed that most of the genes have low or average score, which means that are very ubiquitous.
#' Only a very small fraction of them are season specific (score between 0.9 and 1), around 20 genes.
#'
#'  ## Boxplot gene expression for each gene season specificity score category
par(mar = c(7.1, 7.1, 3.1, 1.1), mgp = c(6, 1, 0))
boxplot(split(es[,2:6],cut(es[,"score"],10)),las=2, main= "Gene average expression relatively to gene season specificity score", 
        ylab = "Gene average expression", xlab = "Gene season specificity score",cex.axis=0.9) 
# "es[,2:6]" is the average expression in individual seasons: autumn, early-summer, late-summer, spring, winter.
# Then it splits genes'average expression in individual seasons relatively to gene's specificity score level organised in 10 bins.
# Finally, it plots genes'average expression for every genes'season specificity score level.

#' Interestingly the genes with less season specificity (ie the most ubiquitous) are the most expressed.
#'
#' ## Check how many genes are never expressed in a season
pander(round(colSums(is.na(es[,2:6])) / nrow(es) * 100,2))
# "colSums(is.na(es[,2:6])" counts all the non expressed genes in a season (their average expression value is listed as na).
# The number of non expressed genes is then divided by the number of genes and multiplied by 100 to calculate the percentage.

#' between 10 and 15% of those genes are never expressed in autumn late summer and winter.
#' Still quite low, not very significative.
#' 
#' ## Check for genes more expressed during a season
#' For every season, select genes whose average expression in a season is equal to the maximal average expression across all seasons.
#' This list of genes obtained corresponds to the genes with an expression peak in the season tested.
#' It also plots gene season specificity distribution per season.

dev.null <- sapply(levels(tp),function(x){
  sel <- es[,x] == es[,"maxn"]
  message(sprintf("There are %s genes with an expression peak in %s",
                  sum(sel,na.rm=TRUE),x))
})

#' Specificity results for cytoscape cluster enriched in photosynthetic genes 2,2:
#' There are 152 genes with an expression peak in early-summer
#' There are 75 genes with an expression peak in late-summer
#' There are 13 genes with an expression peak in autumn
#' There are 59 genes with an expression peak in winter
#' There are 14 genes with an expression peak in spring
#'
#' Specificity results  for cytoscape cluster enriched in photosynthetic genes 2,4:
#' There are 65 genes with an expression peak in early-summer
#' There are 50 genes with an expression peak in late-summer
#' There are 25 genes with an expression peak in autumn
#' There are 17 genes with an expression peak in winter
#' There are 38 genes with an expression peak in spring
#' 
#' Interestingly, cluster 2,4 has more genes specific for Autumn and less for winter than cluster 2,2
#' It could also be interesting to look at the photosynthetic genes that are more expressed in Autumn and Winter.
#'
#' ## Get lists of genes more expressed during a season
season_list <- sapply(levels(tp), function(x){
  sel <- es[,x] == es[,"maxn"]
  sel <- sel[!is.na(sel)]
  gene_list <- names(sel[sel==TRUE])
  return(gene_list)
})

#' Save season specific lists of genes
dir.create(file.path("analysis/season/", "cytoscape_cluster_photosynthetic_genes_enriched_2"))
dev.null <- sapply(levels(tp), function(x){
  write(season_list[[x]],file = paste0("analysis/season/cytoscape_cluster_photosynthetic_genes_enriched_2/",x,"_peak_genes"))
})

#' The lists of season specific genes are too large.
#' There is a need for a threshold to focus only of the most season specific genes with the highest scores.
#'
#' ## Get lists of genes more expressed during a season with high specificity

#' Get the names of genes that are more expressed in a season than in the other, 
#' and whose score is above a threshold specificity score of 0.6.
season_specific_genes <- sapply(levels(tp), function(x){
  x.genes <- names(which(es[es[,x] == es[,"maxn"],"score"]>0.6))
  message(sprintf("There are %s genes that are %s specific",length(x.genes),x))
  return(x.genes)
})

#' High specificity results for cytoscape cluster enriched in photosynthetic genes 2,2:
#' There are 35 genes that are early-summer specific
#' There are 2 genes that are late-summer specific
#' There are 0 genes that are autumn specific
#' There are 5 genes that are winter specific
#' There are 0 genes that are spring specific
#' 
#' High specificity results for cytoscape cluster enriched in photosynthetic genes 2,4:
#' There are 40 genes that are early-summer specific
#' There are 0 genes that are late-summer specific
#' There are 0 genes that are autumn specific
#' There are 0 genes that are winter specific
#' There are 0 genes that are spring specific
#'
#' There are no genes highly specific for Autumn or Winter in cluster 2,4
#' The threshold could be to high for such a short list of genes


#'Example of venn diagram
# ## Overlap
# plot.new()
#grid.draw(venn.diagram(list(Phloem=phloem.winter.genes,
                       #      Xylem=xylem.winter.genes),
                       # filename = NULL,fill=pal[1:2],alpha=0.6))
                       # 
# ### Phloem specific
# write(setdiff(phloem.winter.genes,xylem.winter.genes),file="analysis/season/phloem-only.winter.genes")
# 
# ### Xylem specific
# write(setdiff(xylem.winter.genes,phloem.winter.genes),file="analysis/season/xylem-only.winter.genes")
# 
# ### Common
#write(intersect(phloem.winter.genes,xylem.winter.genes),file="analysis/season/phloem-xylem.winter.genes")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```