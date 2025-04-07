# set wd to source file location
# load the accession list
acc_list <- read.delim("raw_data_names_short.txt", header = FALSE)


mapp_percent_all_samples <- c()

for (i in 1:nrow(acc_list)){
  path <- file.path("quants", acc_list[i,], "logs","salmon_quant.log")
  print(paste0("working on the sample ",  acc_list[i,]))
  log_file <- readLines(path)
  line_with_mapping_rate <- grep("Mapping rate =", log_file)
  mapp_percent <- read.delim(path,
                             header = FALSE, sep = " ",
                             skip = line_with_mapping_rate-1, nrows = 1)[,8]
  mapp_percent_all_samples <- c(mapp_percent_all_samples,
                                as.numeric(substr(mapp_percent,
                                                  start = 1,
                                                  stop = regexpr("%", mapp_percent)-1)))
  rm(mapp_percent, line_with_mapping_rate)
}

combined_summaries <- data.frame(SRR_accession = acc_list$V1,
                                 Mapping_rate = mapp_percent_all_samples)

pdf("histogram_mapping_rate.pdf")
hist(combined_summaries$Mapping_rate,
     main = "Distribution of mapping rates",
     xlim = c(0, 100),
     xlab = "% mapping",
     ylab = "nr. of samples",
     col = "#e5f5e0",
     border = "#31a354")
dev.off()

write.table(combined_summaries, file = "summary_mapping_rate.txt", sep = "\t", row.names = FALSE)
