library(gplots)
library(RColorBrewer)

plot_concordance_heatmap <- function(filename, loh_chromosomes) {
  data <- read.table(filename, sep="\t", header=T, row.names=1)
  colors <- c(rep("white", 22))
  colors[loh_chromosomes] <- "black"
  
  heatmap.2(t(as.matrix(data)), trace="none", scale="none", dendrogram="col", Rowv=F,
            RowSideColors = colors, col=brewer.pal(n = 9, name = "Blues"), 
            key.title="SNP Concordance", 
            key.xlab="% Concordance", 
            key.ylab="", density.info="none", 
            ylab="Chromosome", 
            xlab="Sample pair", srtCol=45, margins=c(8,5))
}

plot_concordance_heatmap_pair <- function(filename, loh_chromosomes) {
  data <- read.table(filename, sep="\t", header=T, row.names=1)
  data[2,] <- data[1,]
  colors <- c(rep("white", 22))
  colors[loh_chromosomes] <- "black"

  heatmap.2(t(as.matrix(data)), trace="none", scale="none", dendrogram="none", Rowv=F,
            RowSideColors = colors, col=brewer.pal(n = 9, name = "Blues"), 
            key.title="SNP Concordance", 
            key.xlab="% Concordance", 
            key.ylab="", density.info="none", 
            ylab="Chromosome", 
            xlab="", srtCol=45, margins=c(8,5), labCol = c("",""))
}


plot_concordance_heatmap("ACC1_matrix.txt", c(1,2,6,10,11,13,14,15,17,18,21,22))
plot_concordance_heatmap_pair("ACC2_matrix.txt", c(2,3,6,8,9,10,11,13,16,18,21,22))
plot_concordance_heatmap_pair("ACC3_matrix.txt", c(1,2,3,6,8,9,11,13,17,18,22))
plot_concordance_heatmap_pair("ACC4_matrix.txt", c(1,2,3,6,8,11,13,15,17,18,22))
plot_concordance_heatmap_pair("ACC5_matrix.txt", c(1,2,3,6,11,12,13,14,15,17,18,22))
plot_concordance_heatmap("ACC6_matrix.txt", c(11, 15, 22)) #+ partial chromosomes
plot_concordance_heatmap_pair("ACC7_matrix.txt", c(1,2,3,4,8,9,10,11,13,15,17,18,21,22))
plot_concordance_heatmap("ACC8_matrix.txt", c(22))
plot_concordance_heatmap("ACC9_matrix.txt", c(1,2,3,8,9,10,11,14,15,17,18,19,20,21,22))
