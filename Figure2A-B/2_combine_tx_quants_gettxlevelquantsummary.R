
# set wd to source file location
# load the accession list
acc_list <- read.delim("raw_data_names_short.txt", header = FALSE)


for (i in 1:nrow(acc_list)){
  path <- file.path("quants", acc_list[i,], "quant.sf")
  quant_results <- read.delim(path, header = TRUE, sep = "\t")
  quant_results <- quant_results[,c(1,4)]
  colnames(quant_results) <- c("Gene", paste("tpm", acc_list[i,], sep = "_"))
  
  if (i==1) {quant_results_all <- quant_results}
  else {quant_results_all <- merge.data.frame(quant_results_all, quant_results, by = "Gene", all = TRUE)}
  rm(quant_results)
}

write.table(quant_results_all, file = "all_quant_tx_Junli.txt", sep = "\t", row.names = FALSE, quote = FALSE)


