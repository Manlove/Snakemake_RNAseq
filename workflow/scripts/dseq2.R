options(repos = c(CRAN = "https://cran.rstudio.com"))
library(DESeq2)
library(dplyr)

run_dseq2 = function(count_data, meta_data, out_path_dseq2, out_path_gsea) {
  data <- read.table(count_data,skip=1,header=TRUE,row.names = 1)
  count_data = data[,c(6:ncol(data))]
  metadata = read.csv(meta_data,header = TRUE)
  
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 1,]
  dds <- DESeq(dds)
  res <- results(dds)
  
  #summary(res)
  #plotMA(res,ylim = c(-2,2))
  #plot(res$log2FoldChange, -log10(res$padj), pch = 20, main = "Volcano Plot",xlab = "Log2 Fold Change", ylab = "-log10 p-value")
  
  write.csv(res,file = out_path_dseq2)
  
  res_df <- as.data.frame(res) %>% filter(!is.na(padj)) %>% mutate(gsea = log2FoldChange * -log10(padj)) %>% select("gsea")
  
  write.table(res_df, out_path_gsea,col.names = TRUE, sep = "\t")
}

run_dseq2(snakemake@input[[1]], snakemake@input[[2]], snakemake@output[[1]], snakemake@output[[2]])
