#' ---
#' title: "Mfuzz Clustering"
#' author: "Thomas Riquelme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```

#' libraries
suppressPackageStartupMessages(library(Mfuzz))

#' load table of counts vst transformed
load(file = "count_data/vst_aware_filtered.rda")

#' # Prepare data

#' transform dates as factor to use to split data for the same date after 
factor_dates <- factor(ordered_dates)

#' merge biological replicates by doing the mean
vst_bio_rep_mean <- do.call(cbind,lapply(split.data.frame(t(vst_aware_filtered),factor_dates),colMeans))

#' Merge biological replicates by doing the median
suppressPackageStartupMessages(library(matrixStats))
vst_bio_rep_median <- do.call(cbind,lapply(split.data.frame(t(vst_aware_filtered),factor_dates),colMedians))

#' test to choose between mean and median
suppressPackageStartupMessages(library("LSD"))
heatscatter(as.vector(vst_bio_rep_mean),as.vector(vst_bio_rep_median))
#' for 99% no difference, the mean seems to have less outliers

#' # Soft clustering (Mfuzz) on photosynthetic gene of interest

#' Aim: To identify expression patterns 

#'read Gene Of Interest csv file
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_list.csv", header = FALSE)
#'read GOI names file
goi_names <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_names.txt", header=FALSE)
#' improve gene of interest data
goi <- cbind(goi,goi_names)
colnames(goi) <- c("id","name")

#' select mean vst counts (of biological replicates for one time point) for gene of interest (goi)
vst_bio_rep_mean_goi <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) %in% goi$id,]

#' Normalize data to obtain z-score to quantify only the variation around the mean (=pattern) and not the amplitude anymore for each each gene
# by doing so we can compare the genes relatively to their expression pattern and not relatively of their amplitude like before
# vst_bio_rep_mean_goi_scaled <- t(scale(t(vst_bio_rep_mean_goi)))
#' ==> not necessary standardisation step will be done with Mfuzz

#' Create the eSet
eset <- ExpressionSet(as.matrix(vst_bio_rep_mean_goi))

#' Standardise
eset.s <- standardise(eset)

#' Estimate the fuzzification
m <- mestimate(eset.s)

#' Find the clusters (12 is based on the previous heatmap)
cl <- mfuzz(eset.s,centers=4,m=m)

mfuzz.plot2(eset.s,
           cl=cl,
           mfrow=c(2,2),
           time.labels = colnames(eset.s),
           las = 2,
           x11=FALSE)

#' get the list of genes for each pattern
acore.list <- acore(eset.s,cl=cl)

temp <- sapply(acore.list,'[[',"NAME")
temp2 <- sapply(temp, as.character)

#' write list of gene from an Mfuzz cluster of interest
write(temp2[[2]],file = "analysis/mfuzz/photosynthetic_goi/winter_peak_listofgenes.txt")

#' # Soft clustering (Mfuzz) on cluster 2,2 enriched for photosynthetic genes 

#' Aim: To compare with photosynthetic goi expression pattern

#'read Gene Of Interest csv file
goi2 <- read.csv(file = "~/Git/UPSCb/projects/spruce-needles/doc/cluster_goi_2,2_gene_list.csv", header = FALSE)

#' select mean vst counts (of biological replicates for one time point) for gene of interest (goi)
vst_bio_rep_mean_goi2 <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) %in% goi2$V1,]

#' Create the eSet
eset2 <- ExpressionSet(as.matrix(vst_bio_rep_mean_goi2))

#' Standardise
eset2.s <- standardise(eset2)

#' Estimate the fuzzification
m2 <- mestimate(eset2.s)

#' Find the clusters (12 is based on the previous heatmap)
cl2 <- mfuzz(eset2.s,centers=4,m=m2)

mfuzz.plot2(eset2.s,
            cl=cl2,
            mfrow=c(2,2),
            time.labels = colnames(eset2.s),
            las = 2,
            x11=FALSE)

#' get the list of genes for each pattern
acore.list <- acore(eset2.s,cl=cl2)

temp <- sapply(acore.list,'[[',"NAME")
temp2 <- sapply(temp, as.character)
#' write list of gene from an Mfuzz cluster of interest
write(temp2[[1]],file = "analysis/mfuzz/cytoscape_cluster_photosynthetic_genes_enriched_1/spiky_cluster1_listofgenes.txt")
write(temp2[[3]],file = "analysis/mfuzz/cytoscape_cluster_photosynthetic_genes_enriched_1/spiky_cluster2_listofgenes.txt")
write(temp2[[2]],file = "analysis/mfuzz/cytoscape_cluster_photosynthetic_genes_enriched_1/not_spiky_cluster1_listofgenes.txt")
write(temp2[[4]],file = "analysis/mfuzz/cytoscape_cluster_photosynthetic_genes_enriched_1/not_spiky_cluster2_listofgenes.txt")

#' # Soft clustering (Mfuzz) on cluster 2,4 enriched with photosynthetic genes 

#' Aim: To compare with photosynthetic goi expression pattern

#'read Gene Of Interest csv file
goi3 <- read.csv(file = "~/Git/UPSCb/projects/spruce-needles/doc/cluster_goi_2,4_gene_list.csv", header = FALSE)

#' select mean vst counts (of biological replicates for one time point) for gene of interest (goi)
vst_bio_rep_mean_goi3 <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) %in% goi3$V1,]

#' Create the eSet
eset3 <- ExpressionSet(as.matrix(vst_bio_rep_mean_goi3))

#' Standardise
eset3.s <- standardise(eset3)

#' Estimate the fuzzification
m3 <- mestimate(eset3.s)

#' Find the clusters (12 is based on the previous heatmap)
cl3 <- mfuzz(eset3.s,centers=4,m=m3)

mfuzz.plot2(eset3.s,
            cl=cl3,
            mfrow=c(2,2),
            time.labels = colnames(eset3.s),
            las = 2,
            x11=FALSE)

#' # Interpretation
#' Cluster 1 of Cytoscape's cluster 2,4 matches cluster 1 of photosynthetic goi.
#' Clusters 2 and 4 of Cytoscape's cluster 2,4 have resemblance in their pattern and they match cluster 4 of photosynthetic goi.
#' And cluster 3 of Cytoscape's cluster 2,4 matches approximately with cluster 3 of photosynthetic goi.
#'
#' In conclusion, we can say that for the two clusters enriched for photosynthetic genes (calculated and visualized in Cytoscape):
#' The first one, cluster 2,2 had most of its genes that match to the expression pattern of cluster 1 in photosynthetic goi, 
#' and also for some of its genes a tendancy to match cluster 2. Plus, there are a new pattern with peaks during summer.
#' In comparison, cluster 2,4 matched less the clusters 1 of goi but matched cluster 3 and 4 of goi.
#'
#' This is logical because in the gene network, the photosynthetic genes were mostly split between these two clusters.
#' We can deduce from that the photosynthetic genes cluster for most of them by expression profile. 
#' Genes whose expression pattern follow pattern 1 and 2 cluster in cluster 2,2 in the network.
#' And genes whose expression pattern follow pattern 3 and 4 cluster in cluster 2,4 in the network.
#' It is interesting to note that in cluster 2,2 there is a other pattern distinct from those of the photosynthetic genes.
#' This other pattern must be linked with the partner genes which have clustered with the photosynthetic genes.
#' The spiky pattern could suggest an activation role. 
#' This could be interesting to look up for those genes and their individual roles to understand why this pattern. 

 
#' ```{r empty,echo=FALSE,eval=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```