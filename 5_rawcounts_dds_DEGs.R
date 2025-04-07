
#dds based on raw counts
# The DESeqDataSet, column metadata, and the design formula
#count matrix called cts and a table of sample information called coldata
library("tximport")
library("DESeq2")
library("gridExtra")
library("ggplot2")
library("ggrepel")
library("edgeR")
library("dplyr")

##DATA####
#read in count data
coldata <- read.delim("JunliRunTable.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = T)
str(coldata)  #printoutput
head(coldata) #printoutput
cts <- as.matrix(read.delim("salmon_counts_counts.txt" , header = TRUE, sep = "\t",row.names = 1))
head(cts)     #printoutput

#create grouped DGEList
#make sure that order of colnames(cts) and coldata$shortname matches

if(all(colnames(cts) == coldata$SRR_txi)){
groups = as.character(coldata$Short_name)
countsDGE = DGEList(round(cts), group = groups)
} else {
  print("colData SSR_txi order is not the same as counttable coldata")
  exit()}


dds <- DESeqDataSetFromMatrix(countData = countsDGE,
                              colData = coldata,
                              design= ~ Genotype + Treatment + Genotype:Treatment)


pdf("parameters_before_filter.pdf", width = 7, height = 5)
hist((log2(counts(dds))), 50,
     xlab = "log2 raw counts",
     main = "Count Distribution before Filtering", col = "#bdd7e7")

hist((colSums(assay(dds)) / 1e6), 30,
     xlab = "M reads per sample",
     main = "Sample Size before Filtering", col = "#bae4b3")
dev.off()
#dds <- estimateSizeFactors(dds)
#counts(dds, normalized=TRUE) #salmon already did the normalization
#dds
#dim(dds)

##FILTERING####
# discard genes with low expression
# 1.) count the amount of samples with expression above count.threshold = 5
# 2.) define the minimum amount of samples (sample.threshold) requiring to pass the threshold1
print(paste("Number of genes before filtering:", nrow(dds)))
count.threshold = 5
sample.threshold = 3
keep <- rowSums(counts(dds) >= count.threshold) >= sample.threshold

print(paste("Percentage of genes with at least",
            count.threshold,"counts in at least", 
            sample.threshold, "samples:",
            round(100*(which(keep %in% TRUE) %>% length()/nrow(dds)), digits = 4),"%"))

dds <- dds[keep,]
print(paste("Number of genes after filtering:", nrow(dds)))

pdf("parameters_after_filter.pdf", width = 7, height = 5)
hist((log2(counts(dds))), 50,
     xlab = "log2 raw counts",
     main = "Count Distribution after Filtering", col = "#6baed6")

hist((colSums(counts(dds)) / 1e6), 30,
     xlab = "M reads per sample",
     main = "Sample Size after Filtering", col = "#74c476")
dev.off()

# data transformation for analysis and visualization
# data for DE analysis stays untransformed

# both vst and rlog functions provide normalization:
# vst recommended for larger datasets (n > 30)，Variance stabilizing transformation (VST) 
# rlog recommended for small datasets (n < 30)
vsd <- vst(dds, blind = TRUE)

pdf("parameters_after_normalization.pdf", width = 7, height = 5)
hist(assay(vsd), 50,
     xlab = "normalized counts",
     main = "VST-normalized Counts", col = "salmon")
dev.off()

# export vst transformed counts
write.table(assay(vsd), "vst_normalized_counts.txt", quote = FALSE)
# color based on different parameters using argument "intgroup"
pca <- plotPCA(vsd, intgroup = c("Short_name"), ntop = nrow(vsd))
pca
# combine plots and export as .pdf
pdf("PCA.pdf", width = 10, height = 10, useDingbats = FALSE)
pca
dev.off()



##REFORMAT GROUP COLUMNS####
dds$group = factor(paste0(dds$Genotype, "_", dds$Treatment))
design(dds) = ~ group

##DeSEQ2####
# calculate differential expression based on reduced design matrix
dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients


##PLOTTING####

#plot comparisons
#define grouping
contrast = c("group", "WT_XopQ", "WT_EV")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame
summary(res)
dim(res)
head(res)

res$gene = row.names(res) #create a gene column
head(res)

#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # up-DEG
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # down-DEG
res_total <- rbind(res_up,res_down)
n_DEG_up <- length(res_up$padj)
n_DEG_up
n_DEG_down <- length(res_down$padj)
n_DEG_down

# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
                  c("log2FoldChange", "padj")]
dim(degs)
head(degs)
# write results_DEs
write.table(degs, "WT_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "WT_UP_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "WT_DOWN_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)


# 依次按照pvalue值log2FoldChange值进行排序
#res <- res[order(res$pvalue, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#plot volcano
v0 <- 
  ggplot(as.data.frame(res), 
         aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_hline(yintercept = -log10(0.05), col = "black",linetype = "dashed") +
  geom_vline(xintercept = -1, col = "black",linetype = "dashed") +
  geom_vline(xintercept = 1, col = "black",linetype = "dashed") +
  geom_point(
    alpha = ifelse(
      res$padj < 0.05 &
        res$log2FoldChange <= -1 |
        res$padj < 0.05 &
        res$log2FoldChange >= 1,
      0,
      0.3),
    shape = 16,
    size = 0.75,
    col = "grey"
  ) +
  #print significant points in red for up and blue for down 
  geom_point (
    alpha = ifelse(
      res$padj < 0.05 &
        res$log2FoldChange <= -1 |
        res$padj < 0.05 &
        res$log2FoldChange >= 1,
      0.5,
      0),
    shape = 16,
    size = 1,
    col = ifelse(
      res$log2FoldChange >= 1,
      "red",
      "blue") 
  )
v0

v1 <- v0 + geom_text_repel(data=head(res, 10), aes(label=gene), 
                           max.overlaps=Inf,
                           min.segment.length = 0) #adding text for the top 10 genes
#ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk
v1


#define grouping
contrast = c("group", "WT_XopQ", "WT_Mock")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame
summary(res)
dim(res)
head(res)

res$gene = row.names(res) #create a gene column
head(res)

#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # up-DEG
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # down-DEG
res_total <- rbind(res_up,res_down)
n_DEG_up <- length(res_up$padj)
n_DEG_up
n_DEG_down <- length(res_down$padj)
n_DEG_down

# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
head(degs)
# write results_DEs
write.table(degs, "WT_XopQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "WT_UP_XopQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "WT_DOWN_XopQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)


# 依次按照pvalue值log2FoldChange值进行排序
#res <- res[order(res$pvalue, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#plot volcano
v0 <- 
  ggplot(as.data.frame(res), 
         aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_hline(yintercept = -log10(0.05), col = "black",linetype = "dashed") +
  geom_vline(xintercept = -1, col = "black",linetype = "dashed") +
  geom_vline(xintercept = 1, col = "black",linetype = "dashed") +
  geom_point(
    alpha = ifelse(
      res$padj < 0.05 &
        res$log2FoldChange <= -1 |
        res$padj < 0.05 &
        res$log2FoldChange >= 1,
      0,
      0.3),
    shape = 16,
    size = 0.75,
    col = "grey"
  ) +
  #print significant points in red for up and blue for down 
  geom_point (
    alpha = ifelse(
      res$padj < 0.05 &
        res$log2FoldChange <= -1 |
        res$padj < 0.05 &
        res$log2FoldChange >= 1,
      0.5,
      0),
    shape = 16,
    size = 1,
    col = ifelse(
      res$log2FoldChange >= 1,
      "red",
      "blue") 
  )
v0

v2 <- v0 + geom_text_repel(data=head(res, 10), aes(label=gene), 
                           max.overlaps=Inf,
                           min.segment.length = 0) #adding text for the top 10 genes
#ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk
v2


#DESeq2_plotCounts show gene expression
data <- plotCounts(dds, gene = "Niben101Scf01212g02011", intgroup = c("Short_name"), returnData = TRUE)   # 指定某个基因, #

#Turn your 'Short_name' column into a character vector
data$Short_name <- as.character(data$Short_name)
#Then turn it back into a factor with the levels in the correct order
data$Short_name <- factor(data$Short_name, levels=unique(data$Short_name))

ggplot(data, aes(x=interaction(Short_name), y=count, color=Short_name)) + 
  geom_jitter(size=2) + 
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("counts (normalized by size factor)")+
  ggtitle("NbGLR2.9a")

data <- plotCounts(dds, gene = "Niben101Scf08670g00020", intgroup = c("Short_name"),returnData = TRUE)   # 指定某个基因
#Turn your 'Short_name' column into a character vector
data$Short_name <- as.character(data$Short_name)
#Then turn it back into a factor with the levels in the correct order
data$Short_name <- factor(data$Short_name, levels=unique(data$Short_name))
ggplot(data, aes(x=interaction(Short_name), y=count, color=Short_name)) + 
  geom_jitter(size=2) + 
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("counts (normalized by size factor)")+
  ggtitle("NbGLR2.9b")

data <- plotCounts(dds, gene = "Niben101Scf00090g03004", intgroup = c("Short_name"),returnData = TRUE)   # 指定某个基因
#Turn your 'Short_name' column into a character vector
data$Short_name <- as.character(data$Short_name)
#Then turn it back into a factor with the levels in the correct order
data$Short_name <- factor(data$Short_name, levels=unique(data$Short_name))
ggplot(data, aes(x=interaction(Short_name), y=count, color=Short_name)) + 
  geom_jitter(size=2) + 
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("counts (normalized by size factor)")+
  ggtitle("NbGLR3.1")


#define grouping
contrast = c("group", "nrg_XopQ", "nrg_EV")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
n_DEG_down <- length(res_down$padj)
n_DEG_down
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "nrg1_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "nrg1_UP_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "nrg1_DOWN_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)


#define grouping
contrast = c("group", "epss_XopQ", "epss_EV")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "epss_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "epss_UP_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "epss_DOWN_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "epssna_XopQ", "epssna_EV")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "epssna_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "epssna_UP_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "epssna_DOWN_XopQ_vs_EV_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)



##################EV_MgCl2
#define grouping
contrast = c("group", "epssna_EV", "epssna_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "epssna_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "epssna_UP_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "epssna_DOWN_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "epss_EV", "epss_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "epss_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "epss_UP_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "epss_DOWN_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "nrg_EV", "nrg_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "nrg1_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "nrg1_UP_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "nrg1_DOWN_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "WT_EV", "WT_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "WT_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "WT_UP_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "WT_DOWN_EV_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)



###################XopQ_Mock
#define grouping
contrast = c("group", "WT_XopQ", "WT_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "WT_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "WT_UP_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "WT_DOWN_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "nrg_XopQ", "nrg_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "nrg1_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "nrg1_UP_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "nrg1_DOWN_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "epss_XopQ", "epss_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # up-DEG
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # down-DEG
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "epss_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "epss_UP_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "epss_DOWN_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

#define grouping
contrast = c("group", "epssna_XopQ", "epssna_MgCl2")
#extract results
res = as.data.frame(results(dds, alpha = 0.05, contrast = contrast))# convert matrix in data frame

res$gene = row.names(res) #create a gene column
head(res)
#create a score combined for log2FoldChange and padj value
# 1.) rank the ABSOLUTE log2fold change
res$rank_l2fc = rank(-abs(res$log2FoldChange))
# 2.) rank the -log(padj) 
res$rank_padj = rank(-abs(-log(res$padj)))
# 3.) combine both scorings and take the highest ranks for each
res$score = apply(res[,c("rank_l2fc", "rank_padj")], 1, FUN = min) 

#order by score
res = res[order(res$score, decreasing = FALSE),]

res_up<- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]      # 表达量显著上升的基因
res_down<- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]    # 表达量显著下降的基因
res_total <- rbind(res_up,res_down)
n_DEG <- length(res_up$padj)
n_DEG
# extract gene names of DEGs
degs <- res[which(abs(res$log2FoldChange) >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)

# write results
write.table(degs, "epssna_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_up-regulated genes
write.table(degs, "epssna_UP_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)

degs <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),
            c("log2FoldChange", "padj")]
dim(degs)
# write results_down-regulated genes
write.table(degs, "epssna_DOWN_XoPQ_vs_Mock_p05_lfc1.txt", sep = "\t",
            quote = FALSE,
            col.names = TRUE, row.names = TRUE)