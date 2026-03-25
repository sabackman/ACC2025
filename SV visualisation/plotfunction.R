library(circlize)
library(data.table)

events_overlap <- function(start1, end1, start2, end2) {
  if (start1 <= end2 && end1 >= start2) {
    return(TRUE)
  }
  FALSE
}
layerEvents <- function(data) {
  data$levels <- rep(1, length(data$Chrom))
  
  for (i in 1:length(data$levels)) {
    j <- 1
    while (j <= length(data$levels)) {
      if (i == j) {
        j <- j+1
        next
      }
      if (data[i,]$Chromosome != data[j,]$Chromosome) {
        j <- j+1
        next
      }
      if (data[i,]$levels != data[j,]$levels) {
        j <- j+1
        next
      }
      if (events_overlap(data[i,]$Start, data[i,]$End, data[j,]$Start, data[j,]$End)) {
        data[i,]$levels <- data[i,]$levels+1
        j <- 1
        next
      }
      j <- j+1
    }
  }
  data
}

plot_genome <- function(vcfpath, titletext="") {
  circos.initializeWithIdeogram(species="hg38")
  
  data <- fread(vcfpath, skip="#CHROM")
  colnames(data)[colnames(data)=="#CHROM"]<-"Chrom"
  
  type <- rep("-1", length(data$Chrom))
  for (i in 1:length(type)) {
    type[i] <- gsub("SVTYPE=", "", strsplit(data[i,]$INFO, ";")[[1]][1])
  }
  data$type <- type
  
  data_dupdel <- data[data$type %in% c("DUP", "DEL"),]
  if (length(data_dupdel$Chrom) > 0) {
    end <- rep(-1, length(data_dupdel$INFO))
    for (i in 1:length(end)) {
      if (!grepl("END=", strsplit(data_dupdel$INFO, ";")[[i]][2])) {
        message("WARNING, format is not as expected. Do not trust results")
      }
      end[i] <- as.integer(gsub("END=", "", strsplit(data_dupdel$INFO, ";")[[i]][2]))
    }
    data_dupdel$end <- end
  
    data_dupdel <- data_dupdel[,c("Chrom", "POS", "end", "type")]
    colnames(data_dupdel) <- c("Chromosome", "Start", "End", "type")
    
    data_dupdel <- layerEvents(data_dupdel)
    
    #Fix coordinates and make sure layer 1 is the outermost one (makes sense for sparsely altered chromosomes)
    data_dupdel$levels <- 1-data_dupdel$levels/max(data_dupdel$levels)
    circos.genomicTrackPlotRegion(
      data_dupdel,
      ylim = c(0, 1),
      track.height = 0.3,
      bg.border = "black",
      stack=F,
      panel.fun = function(region, value, ...) {
        type <- value[[1]]
        layer <- value[[2]]
        circos.genomicLines(
          region,
          layer,
          col = ifelse(type=="DEL", "red", "blue"),
          border = NA,
          type="segment",
          area=F
        )
      }
    )
  }
  
  data_tra <- data[data$type=="TRA",]
  if (length(data_tra$Chrom) > 0) {
    chr2 <- rep("-1", length(data_tra$Chrom))
    end <- rep(-1, length(data_tra$Chrom))
    for (i in 1:length(end)) {
      if (!grepl("CHR2=CHR", strsplit(data_tra[i,]$INFO, ";")[[1]][3])) {
        message("WARNING, format is not as expected. Do not trust results")
      }
      if (!grepl("END=", strsplit(data_tra[i,]$INFO, ";")[[1]][2])) {
        message("WARNING, format is not as expected. Do not trust results")
      }
      
      chr2[i] <- gsub("CHR2=CHR", "chr", strsplit(data_tra[i,]$INFO, ";")[[1]][3])
      end[i] <- as.integer(gsub("END=", "", strsplit(data_tra[i,]$INFO, ";")[[1]][2]))
    }
    data_tra$Chrom2 <- chr2
    data_tra$End <- end
    
    bed1 <- data_tra[,c("Chrom", "POS")]
    bed1$End <- bed1$POS
    bed2 <- data_tra[,c("Chrom2", "End")]
    bed2$POS <- bed2$End
    
    circos.genomicLink(
      bed1,
      bed2)
  }
  circos.clear()
  
  title(titletext)
  
  if (length(data_tra$Chrom)+length(data_dupdel$Chromosome) != length(data$Chrom)) {
    warning("Number of events do NOT add up!!!")
  }
}
