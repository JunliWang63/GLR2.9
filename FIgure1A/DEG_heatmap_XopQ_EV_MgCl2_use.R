setwd("Z:/dep_psl/grp_parker/Junli/WT_nrg1_epss_epssna_RNA_Seq_Junli")
library("pheatmap")
rm(list = ls())
# load the data, sets of DEGs and annotations
vsd <- read.delim("./deseq_output/vst_normalized_counts.txt", header = TRUE, sep = " ")
degs_WT_QE <- row.names(read.table("./deseq_output/WT_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))
degs_nrg1_QE <- row.names(read.table("./deseq_output/nrg1_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))
degs_epss_QE <- row.names(read.table("./deseq_output/epss_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))
degs_epssna_QE <- row.names(read.table("./deseq_output/epssna_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))

degs_WT_QM <- row.names(read.table("./deseq_output/WT_XopQ_vs_Mock_p05_lfc1.txt", header = TRUE))
degs_nrg1_QM <- row.names(read.table("./deseq_output/nrg1_XopQ_vs_Mock_p05_lfc1.txt", header = TRUE))
degs_epss_QM <- row.names(read.table("./deseq_output/epss_XopQ_vs_Mock_p05_lfc1.txt", header = TRUE))
degs_epssna_QM <- row.names(read.table("./deseq_output/epssna_XopQ_vs_Mock_p05_lfc1.txt", header = TRUE))

degs_WT_EM <- row.names(read.table("./deseq_output/WT_EV_vs_Mock_p05_lfc1.txt", header = TRUE))
degs_nrg1_EM <- row.names(read.table("./deseq_output/nrg1_EV_vs_Mock_p05_lfc1.txt", header = TRUE))
degs_epss_EM <- row.names(read.table("./deseq_output/epss_EV_vs_Mock_p05_lfc1.txt", header = TRUE))
degs_epssna_EM <- row.names(read.table("./deseq_output/epssna_EV_vs_Mock_p05_lfc1.txt", header = TRUE))
colData <- read.delim("JunliRunTable.txt", header = TRUE, sep = "\t")

# extract normalized count data for DEGs
degs <- unique(c(degs_WT_EM, degs_WT_QE,degs_WT_QM,degs_nrg1_EM,degs_nrg1_QE, degs_nrg1_QM,
                 degs_epss_EM,degs_epss_QE,degs_epss_QM,
                 degs_epssna_EM, degs_epssna_QE,degs_epssna_QM))
degCounts <- vsd[degs,]
head(rownames(degCounts), 5)
tail(rownames(degCounts), 5)

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
pheatmap(degZScores, cluster_rows=TRUE, show_rownames=FALSE,
         gaps_col=c(4,8), cluster_cols=FALSE, show_colnames = TRUE, cutree_rows = 4,
         border_color = NA, annotation_col=group,
         filename = "heatmap_DEGs_XopQ_EV_MgCl2_use.pdf")
