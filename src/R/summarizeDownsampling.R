
# Correlations of the 3 downsampled datasets

downsampling = c("0.125", "0.25", "0.50", "0.75")
downsampling = c("0.125", "0.25", "0.50")

data.l = list()

for (downSamplingCur in downsampling) {
  
  file = paste0("/scratch/carnold/CLL/TF_act_downsampling", downSamplingCur, "/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.tsv")
  assertFileExists(file)
  
  data.l[[downSamplingCur]] = read_tsv(file, col_names = TRUE)
  
}

file = paste0("/scratch/carnold/CLL/TF_act/output/FINAL_OUTPUT/extension100/CLL.all.volcano.l2FC.removedCGbias.AR.FDR5.tsv")
assertFileExists(file)
data.l[["1"]] = read_tsv(file, col_names = TRUE)

names(data.l) =  c(paste0("fraction", downsampling), "full")

plots.l = list()


for (downSamplingCur in downsampling) {
  
  
  for (typeCur in c("weighted_mean", "weighted_Tstat", "weighted_CD", "weighted_median")) {
    
    
    data.df = as_tibble(map(data.l, typeCur))
    
    
    p <- ggplot(data.df, aes_string(x = paste0("fraction", downSamplingCur), y = "full"))
    p <- p + geom_point()
    p <- p + .getThemeForGGPlot()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p <- p + ggtitle(paste0("Pearson correlation ", typeCur, " ", downSamplingCur, " vs. all: ", round(cor(data.df[, paste0("fraction", downSamplingCur)], data.df$full), 2)))
    p <- p + geom_smooth(method = "loess", size = 1.5)
    plots.l[[paste0(downSamplingCur, "_", typeCur)]] = p
  }
  
}

printMultipleGraphsPerPage(plots.l, nCol = 2, nRow = 2, pdfFile = "downsamplingResults.pdf", width = 12, height = 12) 

