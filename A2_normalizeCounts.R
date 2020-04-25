#' ---
#' title: "Spruce Needles Kallisto Biological QA"
#' author: "Nicolas Delhomme & Thomas Riquelme"
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

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#'#### Get the count data

#' # Raw data
#' ## Loading
#' Read the sample information
#' Create a table from the sequencing data (filename,sequencing date)
samples <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/rna_seq_spruce_needles_2011_2012_sample_info.csv")


#'List all the input kallisto tsv files and store them in orig variable
orig <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(orig) <- sub("_sortmerna.*","",
                    sapply(strsplit(orig, "/"), .subset, 2))
# sub replaces the pattern in arg1 with arg2 in all the elements of vectors in arg3                   

#' Match samples data to samples infos

orig <- orig[match(samples$Sequencing_id,names(orig))]
orig <- orig[!is.na(names(orig))]

samples <- samples[samples$Sequencing_id %in% names(orig),]

stopifnot(all(samples$Sequencing_id == names(orig)))

#' Extract the expression data from kallisto output files and store it into matrices
tx <- suppressMessages(tximport(files = orig, 
                                type = "kallisto", 
                                txOut = TRUE))
# tximport(): looks inside the kallisto output files abundance.tsv, extract information and put it nicely into matrices for further use in downstream analysis
#store the matrices into tx
# str(tx)= a list of 4: abundance matrix, counts matrix, length matrix, countsFromAbundance (chr "no")

#' Store the counts in kt variable
kt <- round(tx$counts)
# round (keep only integer) and store the counts in kt 
# str(kt): a matrix
# rows: 66000 genes
# columns: 166 samples

#' Check for the genes that are never expressed
sel <- rowSums(kt) == 0 
# str(sel)= logical vector
# sums the counts for one gene in each samples
# show TRUE or FALSE for each gene if the gene is never expressed (has 0 counts in every time samples)

#' Check for the number of genes never expressed
sprintf("%s%% (%s) of %s genes/transcripts are not expressed",
        round(sum(sel) * 100/ nrow(kt),digits=1),
        sum(sel),
        nrow(kt))
#in our experiment the number of genes and transcript are the same (one transcript for one gene)

#' Display the samples mean raw counts distribution
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#'
plot(density(log10(rowMeans(kt))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")
#density=likelihood, chance
# plot the likelihood that the mean of genes has this mean of raw counts 
# give an idea of the read depth, here the majority of gene maps ~30 reads (1.5 log10) 

#' Display all samples raw counts distribution
cols <- sample(pal,30,replace = TRUE)
plot.multidensity(mclapply(1:ncol(kt),function(k){log10(kt)[,k]},mc.cores=16L),
                  col=cols[as.integer(samples$Sampling_Date)],
                  legend.x="topright",
                  legend=levels(samples$Sampling_Date),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")
# plot the likelihood that a gene in different time points has a certain number of counts 
# genes have a lot of chances (0.2 to 0.4 depending on the dates) to have a count of one 
# genes also have a lot of chances (0.2 to 0.3) to have between 2 and 100 counts 
# genes have fewer chances (0.2 to 0.0) to have between 100 and a 1000 counts per gene
# genes have even fewer chances (0.05 to 0.0) to have very high counts (1000 to 100000)= "powertail" 

#==> one sample seems to have a lower sequencing depth (lower total of counts compare to the others)

#' find the bad sample with the minimum total of counts per sample 
kt_col_sum <- colSums(kt)
min_kt <- min(kt_col_sum)
bad_sample <- names(kt_col_sum)[kt_col_sum == min_kt]
#' remove from kt 
kt <- kt[,!colnames(kt) %in% bad_sample]
#' and remove from samples
bad_sample2 <- rownames(samples)[samples$Sequencing_id == bad_sample]
samples <- samples[!rownames(samples) %in% bad_sample2,]

#' plot it again
#' plot cumulative gene coverage
plot(density(log10(rowMeans(kt))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")
#' plot gene coverage for each samples
plot.multidensity(mclapply(1:ncol(kt),function(k){log10(kt)[,k]},mc.cores=16L),
                  col=cols[as.integer(samples$Sampling_Date)],
                  legend.x="topright",
                  legend=levels(samples$Sampling_Date),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("analysis","kallisto"),showWarnings=FALSE,recursive=TRUE)
write.csv(kt,file="analysis/kallisto/spruce-needles-unormalised-gene-expression_data.csv")
save(kt, samples, file = "counts.rda")

#load("counts.rda")

#'#### Merge technical replicates

#' Merge technical replicates 

kt_tech_merged <- do.call(cbind,lapply(split.data.frame(t(kt),samples$sampling_id),colSums))
# transpose the count table kt to get the samples' names as rows in order to process them next (functions process rows and not columns)
# then, do the sums of counts for all genes at one date
# finally, bind again the splitted rows to reunified a count table with technical replicates merged for each dates 

#' Merge sampling dates for technical replicates
sampling_dates_merged <- sapply(split(as.character(samples$Sampling_Date),samples$sampling_id),unique)

#' Save count table with technical replicates merged
save(sampling_dates_merged, kt_tech_merged, file="kt_tech_merged.rda")


#'#### Discard the NE samples which were resequenced with higher depth as TR00 (keep TR00 with the same number)

#' find the bad samples
samples_names <- colnames(kt_tech_merged)
samples_names
new_names <- sub("NE-|TR00","",samples_names)
new_names
list <- split(samples_names, new_names)
list
duplicates <- list[sapply(list,length)==2]
bad_samples <- sapply(X = 1:length(duplicates),function(i){duplicates[[i]][1]})

#' remove from kt_tech_merged 
kt_tech_merged_cropped <- kt_tech_merged[,!colnames(kt_tech_merged) %in% bad_samples]

#' and remove the dates for this samples in the dates vector
sampling_dates_merged_cropped <- sampling_dates_merged[!names(sampling_dates_merged) %in% bad_samples]

# #' Save cropped dates and counts 
# Save(kt_tech_merged_cropped, sampling_dates_merged_cropped, file="kt_tech_merged_cropped.rda")

#' order dates chronologically
dates <- sampling_dates_merged_cropped
ordered_dates <- dates[order(as.Date(dates, format="%Y-%m-%d"))]

#' order count table chronologically
kt_tech_merged_cropped_ordered <- kt_tech_merged_cropped[,names(ordered_dates)]

#' Save ordered dates and counts
save(ordered_dates, kt_tech_merged_cropped_ordered, file = "kt_tech_merged_cropped_ordered.rda")

#'##### Data normalisation 

#' ### vst blind

#load technical replicates merged data cropped for NE bad samples
load("kt_tech_merged_cropped_ordered.rda")

#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
#'
#' Create the dds object, without giving any prior on the design

dds.kt <- DESeqDataSetFromMatrix(
  countData = kt_tech_merged_cropped_ordered,
  colData = data.frame(date=ordered_dates),
  design = ~ date)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kt <- estimateSizeFactors(dds.kt)
sizes.kt <- sizeFactors(dds.kt)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kt)
boxplot(sizes.kt, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation

vsd.kt <- varianceStabilizingTransformation(dds.kt, blind=TRUE)
vst_blind <- assay(vsd.kt)
#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kt+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kt,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5

write.csv(vst_blind,"analysis/kallisto/spruce-needles-normalized-data-vst-blind.csv")

save(ordered_dates, vst_blind, file="vst_blind.rda")

#' ### vst aware (calculated accordingly to samplingDates)

#' Load count table with technical replicates merged and cropped for NE bad samples
load("kt_tech_merged_cropped_ordered.rda")
 
#' Calculate vst aware with DESeq2 taking the date as a parameter
#' Create the dds object
dds2 <- DESeqDataSetFromMatrix(
  countData = kt_tech_merged_cropped_ordered,
  colData = data.frame(date = ordered_dates),
  design = ~ date)

#' Check the size factors (i.e. the sequencing library size effect)
dds2 <- estimateSizeFactors(dds2)
sizes2 <- sizeFactors(dds2)
pander(sizes2)
boxplot(sizes2, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation

vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)
vst_aware <- assay(vsd2)
#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kt+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds2,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind

write.csv(vst_aware,"analysis/kallisto/spruce-needles-normalised-data-vst-aware.csv")
save(ordered_dates, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' ### Filter data for expressed enough genes 

#' Select genes that are expressed above a threshold of counts 'exp' and in at least 'n' replicates
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

#' find the right threshold to get around 20000 genes
sapply(X = 1:10, function(i){sum(geneSelect(vst_aware,ordered_dates,i))})

#' selection of genes whose counts are over or equal to a vst value of 1.29 in at least 2 replicates of one time point
catalog.sel <- geneSelect(vst_aware,ordered_dates,1.29)
sum(catalog.sel)

vst_aware_filtered <- vst_aware[catalog.sel,]

save(vst_aware_filtered, ordered_dates, file = "vst_aware_filtered.rda")

write.table(vst_aware_filtered, file=file.path("vala/data","spruce_needles_rnaseq_vst_aware_filtered.tsv"), quote=FALSE, sep="\t")


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```       