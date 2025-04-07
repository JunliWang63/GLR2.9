
library("pheatmap")

# load the data, sets of genes and annotations
vsd <- read.delim("vst_normalized_counts.txt", header = TRUE, sep = " ")
nlr <- read.delim("NbNRG1_NRG2_ESP_geneID.txt", header = TRUE, sep = "\t")
colData <- read.delim("JunliRunTable.txt", header = TRUE, sep = "\t")
head(colData)
head(nlr)
tail(nlr)
# extract normalized count data for genes of interest
degCounts <- vsd[nlr$GeneID,]
degCounts
dim(degCounts)

# make identifiers for combining biological replicates
annotation <- data.frame(SRR = as.factor(colData$Run),
                         ident = as.factor(colData$Short_name))

# combine all biological replicates
templist <- list()

for(i in 1:length(levels(annotation$ident))) {
  
  # id is current Genotype_Timepoint_Treatment combination
  id <- levels(annotation$ident)[i]
  print(paste("Working on:", id, sep = "  "))
  
  SRRs <- annotation$SRR[which(annotation$ident == id)]
  counts <- degCounts[,which(colnames(vsd) %in% SRRs)]
  
  means <- as.data.frame(apply(counts, 1, mean))
  colnames(means) <- id
  
  templist[[i]] <- means
  rm(id, SRRs, counts, means)
}

vsdMeans <- do.call(cbind, templist)
rm(templist)


# perform z-score normalization using scale function
degZScores <- as.data.frame(t(scale(t(vsdMeans), center = TRUE, scale = TRUE)))

# reorder columns in degZScores
table(colData$Short_name)
ordered <- data.frame(Treatment=c("WT_MgCl2", "nrg_MgCl2", "epss_MgCl2", "epssna_MgCl2",
                                  "WT_EV", "nrg_EV", "epss_EV", "epssna_EV",
                                  "WT_XopQ", "nrg_XopQ", "epss_XopQ", "epssna_XopQ"))

col_order <- ordered$Treatment

degZScores <- degZScores[, col_order]


# build the heatmap and make annotation
group=data.frame(type=c(rep("Mock",4),rep("Pf0_EV",4),rep("Pf0_XopQ",4)))

rownames(group)=colnames(degZScores)
pheatmap(degZScores, cluster_rows=TRUE, show_rownames=TRUE,
         gaps_col=c(4,8), cluster_cols=FALSE, show_colnames = TRUE, cutree_rows = 1,
         border_color = NA, annotation_col=group,
         filename = "heatmap_NRG1_2_ESP_DEGs_XopQ_EV_MgCl2.pdf")





