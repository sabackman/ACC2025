library(MutationalPatterns)
library(gplots)
library("BSgenome.Hsapiens.UCSC.hg38", character.only=TRUE)

ref_genome<- "BSgenome.Hsapiens.UCSC.hg38"


samplenames <- c("ACC1_7504", "ACC1_7510", "ACC1_8745", "ACC1_8748", "ACC1_8749", 
                 "ACC2_10817-2", "ACC2_10820-2",
                 "ACC3_5589", "ACC3_6031", 
                 "ACC4_2269", "ACC4_4465", 
                 "ACC5_18381-nr3", "ACC5_11529",
                 "ACC6_11259", "ACC6_11260", "ACC6_8553-nr4", "ACC6_A7551-20AT", 
                 "ACC7_11599", "ACC7_2177-nr5", 
                 "ACC8_11691", "ACC8_11694", "ACC8_11712", "ACC8_A9008-19", "ACC8_11912", "ACC8_11695",
                 "ACC9_11913", "ACC9_11913n2", "ACC9_11914", "ACC9_11504")
mapping <- c("ACC1-P1", "ACC1-P2", "ACC1-M1", "ACC1-M2", "ACC1-M3", 
             "ACC2-P", "ACC2-M", 
             "ACC3-P", "ACC3-R", 
             "ACC4-P", "ACC4-R", 
             "ACC5-P", "ACC5-M",
             "ACC6-P1", "ACC6-P2", "ACC6-M1", "ACC6-M2", 
             "ACC7-P", "ACC7-M", 
             "ACC8-P1", "ACC8-P2", "ACC8-M1", "ACC8-M2", "ACC8-M3", "ACC8-Thrombus",
             "ACC9-P1", "ACC9-P2", "ACC9-R", "ACC9-M")

samplemapping <- data.frame(name=samplenames, mapping=mapping)

vcf_files <- list.files("filteredvcfs", full.names = TRUE, pattern=".vcf")
sample_names <- gsub(".vcf", "", gsub("filteredvcfs/", "", vcf_files))

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)

type_occurrences <- mut_type_occurrences(grl, ref_genome)

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)


#
# Strict refitting
#
signatures <- get_known_signatures()
fit_res_ <- fit_to_signatures_strict(mut_mat, signatures)

fit_res <- fit_res_$fit_res

contributions <- t(fit_res$contribution)


#Rename samples for publication ready figure
newnames <- rownames(contributions)
for (i in 1:length(newnames)) {
  newnames[i] <- samplemapping[samplemapping$name==newnames[i],]$mapping
}
rownames(contributions) <- newnames



normalized <- contributions/rowSums(contributions)

heatmap.2(as.matrix(normalized[,colSums(normalized>=0.1)>0]), trace="none", Rowv=FALSE, dendrogram="col", 
          scale="none", col="bluered",
          density.info="none", key.title=NA, margins=c(5,8), key.xlab="Contribution (fraction)")

plot_original_vs_reconstructed(mut_mat, fitres_final$reconstructed, 
                               y_intercept = 0.95)

#
# Rainfall plots
#
chromosomes <- seqnames(get(ref_genome))[1:22]

# Make a rainfall plot
for (i in 1:29) {
  png(filename=file.path("RainfallPlots", paste(names(grl[i]), ".png", sep="")), width=800, height=500)
  p <- plot_rainfall(grl[[i]],
                title = names(grl[i]),
                chromosomes = chromosomes, cex = 0.7, ylim = 1e+09
  )
  print(p)
  dev.off()
}


