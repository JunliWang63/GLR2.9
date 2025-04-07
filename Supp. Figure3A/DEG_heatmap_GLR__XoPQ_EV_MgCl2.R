
library("pheatmap")

# load the data, sets of genes and annotations
vsd <- read.delim("vst_normalized_counts.txt", header = TRUE, sep = " ")
nlr <- read.delim("GLR_Junli.txt", header = TRUE, sep = "\t")
colData <- read.delim("JunliRunTable.txt", header = TRUE, sep = "\t")
head(nlr)
tail(nlr)
# extract normalized count data for genes of interest
degCounts <- vsd[nlr$GeneID,]
head(rownames(degCounts), 5)
tail(rownames(degCounts), 5)

degCounts
head(degCounts)
tail(degCounts)
sum(is.na(degCounts))
giveNAs = which(is.na(degCounts), arr.ind=TRUE)
head(giveNAs)
#give examples of rows with NA
degCounts[c(6,11,16),]
#We get the NA rows out and start checking what to remove:
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(degCounts[-as.numeric(names(tab[1:i])),]))
})
rmv = names(tab)[1:min(which(checkNA==0))]

degCounts1 = degCounts[-as.numeric(rmv),]
dim(degCounts1)#remove 6 rows with NA

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
  counts <- degCounts1[,which(colnames(vsd) %in% SRRs)]
  
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
         filename = "heatmap_GLR__XopQ_EV_MgCl2.pdf")



