options(repos = c(CRAN = "https://cran.rstudio.com"))

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('PCAtools')
#BiocManager::install('DESeq2')

library(DESeq2)
library(PCAtools)
library(ggplot2)
library(dplyr)


create_pca = function(count_data, meta_data, out_path) {
  counts = read.csv(count_data,sep="\t", skip = 1, header=TRUE, row.names = 1) %>%
    select((6:ncol(.)))
  metadata = read.csv(meta_data, header=TRUE)

  dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~condition)
  dds = DESeq(dds)
  
  vsd = vst(dds, blind = TRUE)
  pca_data = assay(dds)
  p = pca(pca_data, metadata = colData(dds), removeVar = 0.8)
  screeplot(p)
  outplot = biplot(p,
                  colby = "condition",
                  legendPosition = "right",
                  lab = rownames(metadata))
  
  # pca = prcomp(t(assay(dds)))
  # pca_data = as.data.frame(pca$x)
  # pca_data$condition = metadata$condition
  
  # outplot = ggplot(pca_data, aes(x=PC1,y=PC2,color = condition)) + 
  #   geom_point(size = 3) + 
  #   labs(title = "PCA", x = "PC1", y = "PC2") + 
  #   theme_classic() + 
  #   theme(legend.title = element_blank())
  # outplot
  ggsave(filename = out_path, width = 6, height = 4, device = "png")
}

create_pca(snakemake@input[[1]], snakemake@input[[2]], snakemake@output[[1]])