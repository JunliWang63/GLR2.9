install.packages("ggvenn")
suppressPackageStartupMessages(library(ggvenn))
install.packages("VennDiagram")
library(VennDiagram)
# load the data, sets of genes and annotations

degs_WT <- row.names(read.table("WT_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))
degs_nrg1 <- row.names(read.table("nrg1_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))
degs_epss <- row.names(read.table("epss_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))
degs_epssna <- row.names(read.table("epssna_XopQ_vs_EV_p05_lfc1.txt", header = TRUE))

degs_WT1 <- as.data.frame(degs_WT)
rownames(degs_WT1) <- rownames(degs_WT)

degs_nrg11 <- as.data.frame(degs_nrg1)
rownames(degs_nrg11) <- rownames(degs_nrg1)

degs_epss1 <- as.data.frame(degs_epss)
rownames(degs_epss1) <- rownames(degs_epss)

degs_epssna1 <- as.data.frame(degs_epssna)
rownames(degs_epssna1) <- rownames(degs_epssna)
  
a <- list(WT=rownames(degs_WT1),
            nrg1=rownames(degs_nrg11),
            epss1=rownames(degs_epss1),
            epssna1=rownames(degs_epssna1))
ggvenn(a,show_percentage = FALSE, digits = 1,fill_color = c("blue", "yellow", "green", "red"),text_size = 4,
       stroke_alpha = 1,
       stroke_size = 1)

