#' ---
#' title: "Plot Expression Profile with Confidence Intervals"
#' author: "Thomas Riquelme"
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

#Libraries
library(ggplot2)
library(matrixStats)

#load table of counts vst transformed
load(file = "vst_aware.rda")

dates <- ordered_dates
unique_dates <- unique(dates)
factor_dates <- factor(ordered_dates)


vst_bio_rep_mean <- do.call(cbind,lapply(split.data.frame(t(vst_aware),factor_dates),colMeans))

vst_bio_rep_median <- do.call(cbind,lapply(split.data.frame(t(vst_aware),factor_dates),colMedians))
rownames(vst_bio_rep_median) <- rownames(vst_bio_rep_mean)

vst_bio_rep_sd <- do.call(cbind,lapply(split.data.frame(t(vst_aware),factor_dates),colSds))

vst_bio_rep_lwr <- vst_bio_rep_mean-qt(0.975,3)*vst_bio_rep_sd/sqrt(3)

# "2011-06-30" has no replicates, thus it is impossible to calculate sd, lwr and upr
# replacement of NA values  in upr and lwr by the mean to be able to plot the confidence interval of the mean later
vst_bio_rep_lwr[,"2011-06-30"] <- vst_bio_rep_mean[,"2011-06-30"]
vst_bio_rep_lwr[vst_bio_rep_lwr < 0] <- 0

vst_bio_rep_upr <- vst_bio_rep_mean+qt(0.975,3)*vst_bio_rep_sd/sqrt(3)

vst_bio_rep_upr[,"2011-06-30"] <- vst_bio_rep_mean[,"2011-06-30"]


#read Gene Of Interest csv file
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_list.csv", header = FALSE)

#read GOI names file
goi_names <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_names.txt", header=FALSE)

# Test to plot one gene of the list
psbs_mean <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) == goi[52,1]]
psbs_median <- vst_bio_rep_median[rownames(vst_bio_rep_median) == goi[1,1]]
psbs_lwr <- vst_bio_rep_lwr[rownames(vst_bio_rep_lwr) == goi[52,1]]
psbs_upr <- vst_bio_rep_upr[rownames(vst_bio_rep_upr) == goi[52,1]]

DF <- data.frame(time=unique_dates,
                 mean=psbs_mean,
                 lwr=psbs_lwr,
                 upr=psbs_upr)

# ggplot(DF, aes(Time)) + 
#   geom_line(aes(y=menle), colour="blue") + 
#   geom_ribbon(aes(ymin=menlelb, ymax=menleub), alpha=0.2)

p <- ggplot(DF, aes(time, group = 1)) +
  geom_line(aes(y=mean), color="blue") +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "vst counts", title = "New plot title") +
  ggtitle("Expression of gene test", subtitle = NULL)
plot(p)


# do a loop to plot the expression profile for each gene of the list
#n=0
for (i in 1:nrow(goi)) {
  message(i)
  #check if goi are in vst with %in%
  if (goi[i,1] %in% rownames(vst_aware)) {
    g_median <- vst_bio_rep_median[rownames(vst_bio_rep_median) == goi[i,1]]
    g_mean <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) == goi[i,1]]
    g_lwr <- vst_bio_rep_lwr[rownames(vst_bio_rep_lwr) == goi[i,1]]
    g_upr <- vst_bio_rep_upr[rownames(vst_bio_rep_upr) == goi[i,1]]
    DF <- data.frame(time=unique_dates,
                     mean=g_median,
                     lwr=g_lwr,
                     upr=g_upr)
    p <-ggplot(DF, aes(time, group = 1)) +
      geom_line(aes(y=mean), color="blue") +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(y = "vst counts", title = paste("Expression profile of", goi_names[i,1]))
    plot(p)
    ggsave(filename = paste0("expression_profiles/",goi_names[i,1],"_expression_profile.jpeg"))
    #n=n+1
  }
}
#print(n)

#ELIP_C does not match our data genes id

#'### Hierarchical clustering to see which genes have similar patterns

# improve gene of interest data
goi <- cbind(goi,goi_names)
colnames(goi) <- c("id","name")


# select mean vst counts (of biological replicates for one time point) for gene of interest (goi)
vst_bio_rep_mean_goi <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) %in% goi$id,]

# replace gene ids by gene names as rownames
goi <- goi[! goi$name == "ELIP_C",]
goi <- goi[order(match(goi$id,rownames(vst_bio_rep_mean_goi))),]
rownames(vst_bio_rep_mean_goi) <- goi$name

# Normalize data to obtain z-score to quantify only the variation around the mean (=pattern) and not the amplitude anymore for each gene
# by doing so we can compare the genes relatively to their expression pattern and not relatively of their amplitude like before
vst_bio_rep_mean_goi_scaled <- t(scale(t(vst_bio_rep_mean_goi)))
#head(vst_bio_rep_mean_goi_scaled)

# qqplot to check if our genes expression follows a Normal distribution or not
require(graphics)
qqnorm(vst_bio_rep_mean_goi)
qqline(vst_bio_rep_mean_goi, col=3)
qqnorm(vst_bio_rep_mean_goi_scaled)
qqline(vst_bio_rep_mean_goi_scaled, col = 2)
# ==> does not follow a line, the distribution is not Normal
# ==> then Pearson correlation should be use because is not parametric

#Perform hierarchical clustering
#Firstly compute distance according to "spearman correlation" method
#Secondly compute clustering according to "complete" linkage
library(amap)
vst_bio_rep_mean_goi_cluster <- hcluster(vst_bio_rep_mean_goi, method = "spearman", link = "complete")
vst_bio_rep_mean_goi_scaled_cluster <- hcluster(vst_bio_rep_mean_goi_scaled, method = "spearman", link = "complete")

#' Plot the dendrogram
plot(vst_bio_rep_mean_goi_cluster, main="Cluster Dendrogram of photosynthetic genes")
plot(vst_bio_rep_mean_goi_scaled_cluster, main="Cluster Dendrogram of photosynthetic genes")

library(corrplot)

# mar <- par("mar")
# par(mar=c(6.1,4.1,6.1,2.1))

pdf(file = "expression_profiles/photosynthetic_genes_correlation_matrix.pdf")
corrplot(cor(t(vst_bio_rep_mean_goi)),method="square",order="hclust", title = "Correlation matrix of photosynthetic genes",
         tl.cex=0.4, cl.cex=0.5, tl.col="black", addrect=2, is.corr = FALSE, mar = c(0,0,2,0))
dev.off()

# par(mar=mar)

corrplot(cor(t(vst_bio_rep_mean_goi_scaled)),method="square",order="hclust",
        tl.cex=0.4, cl.cex=0.5, tl.col="black", addrect=2, is.corr = FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```