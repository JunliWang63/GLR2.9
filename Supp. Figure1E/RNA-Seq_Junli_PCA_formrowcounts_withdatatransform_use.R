library("vegan")
library("ggplot2")
fname <- file.choose()#choose vst_normalized_counts
fname
pdata <- read.table(fname)
pdata
dim(pdata)
head(pdata)
pdata1 <- apply(pdata, 2, function(x) as.numeric(x))#character_numeric
pdata1 <- as.data.frame(pdata1)
head(pdata1)
rownames(pdata1) <- rownames(pdata)
head(pdata1)

pca <- rda(pdata1, scale = TRUE)
pca
pdatapca <- as.data.frame(pca$CA$v)
head(pdatapca)
coldata <- read.delim("JunliRunTable.txt", header = TRUE, sep = "\t")
str(coldata)
head(coldata)
pdatapca$genotype <- coldata[match(rownames(pdatapca), coldata[,1]), 2]
pdatapca$treatment <- coldata[match(rownames(pdatapca), coldata[,1]), 3]
ggplot(pdatapca, aes(PC1, PC2, color=treatment, shape=genotype))+
  geom_point() +
  ggtitle("PCA")

ggsave("all_samples_PCA.pdf", width = 4, height = 4, useDingbats = FALSE)

  
  
  
      


