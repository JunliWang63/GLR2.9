# preparation of tx2gene for tximport to summarize
# transcript level quantification to the gene level

# BiocManager::install("Biostrings")
library("Biostrings")

# read reference transcript file
ref <- readDNAStringSet(file.path("Niben101_annotation.transcripts.fa"))
ref[1:5]

### extract transcript names
tx_name <- ref@ranges@NAMES

### extract gene name
# in some cases transcriptID is the same as geneID, without isoform info
gene_name <- tx_name
# get transcript indexes that have isoform information
tx_with_isoform_info <- grep("\\.", tx_name)
# extract geneID from transcriptID when the isoform info is available
gene_name[tx_with_isoform_info] <- substr(tx_name[tx_with_isoform_info],
                                          start = 1, stop = regexpr("\\.", tx_name)-1)

### prepare a df tx2gene with columns transcriptID and geneID
tx2gene <- data.frame(transcriptID = tx_name, geneID = gene_name)

### write the tx2gene df into a file
write.table(tx2gene, file = "Z:/dep_psl/grp_parker/Junli/WT_nrg1_epss_epssna_RNA_Seq_Junli/reference/Niben101_annotation.transcripts.fa/tx2gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)
fname <- file.choose()
file.exists(fname)
tx2gene <- read.delim(fname,stringsAsFactors = FALSE)
summary(tx2gene)
dim(tx2gene)
head(tx2gene)
