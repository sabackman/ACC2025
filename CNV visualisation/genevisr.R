library(GenVisR)
library(data.table)

#
# Custom genome without chrX and chrY as we do not provide copy number estimates
# for these
#
custom <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == "hg38", ]
custom <- custom[!custom$chrom %in% c("chrX", "chrY"), c("chrom", "chromStart", "chromEnd")]
colnames(custom) <- c("chromosome", "start", "end")
#ACC1
data1 <- read.table("SJ-2335-7504_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("SJ-2335-8745-v1_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data3 <- read.table("SJ-2335-8748_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data4 <- read.table("SJ-2335-8749-v1_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data5 <- read.table("UI-3077-7510-v1_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC1 <- as.data.table(rbind(data1, data2, data3, data4, data5))
ACC1[ACC1$sample=="SJ-2335-7504",]$sample <- "ACC1-P1"
ACC1[ACC1$sample=="UI-3077-7510-v1",]$sample <- "ACC1-P2"
ACC1[ACC1$sample=="SJ-2335-8745-v1",]$sample <- "ACC1-M1"
ACC1[ACC1$sample=="SJ-2335-8748",]$sample <- "ACC1-M2"
ACC1[ACC1$sample=="SJ-2335-8749-v1",]$sample <- "ACC1-M3"


cnSpec(ACC1, y=custom)

#ACC2
data1 <- read.table("SJ-2335-10817-2_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("SJ-2335-10820-2_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC2 <- as.data.table(rbind(data1, data2))
ACC2[ACC2$sample=="SJ-2335-10817-2",]$sample <- "ACC2-P"
ACC2[ACC2$sample=="SJ-2335-10820-2",]$sample <- "ACC2-M"
cnSpec(ACC2, y=custom)

#ACC3
data1 <- read.table("SJ-2335-5589_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("SJ-2335-6031_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC3 <- as.data.table(rbind(data1, data2))
ACC3[ACC3$sample=="SJ-2335-5589",]$sample <- "ACC3-P"
ACC3[ACC3$sample=="SJ-2335-6031",]$sample <- "ACC3-R"
cnSpec(ACC3, y=custom)

#ACC4
data1 <- read.table("SJ-2335-2269_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("SJ-2335-4465_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC4 <- as.data.table(rbind(data1, data2))
ACC4[ACC4$sample=="SJ-2335-2269",]$sample <- "ACC4-P"
ACC4[ACC4$sample=="SJ-2335-4465",]$sample <- "ACC4-R"
cnSpec(ACC4, y=custom)

#ACC5
data1 <- read.table("SJ-2335-18381-nr3_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("SJ-2335-11529_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC5 <- as.data.table(rbind(data1, data2))
ACC5[ACC5$sample=="SJ-2335-18381-nr3",]$sample <- "ACC5-P"
ACC5[ACC5$sample=="SJ-2335-11529",]$sample <- "ACC5-M"
cnSpec(ACC5, y=custom)

#ACC6
data1 <- read.table("SJ-2335-11259-v1_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("UI-3077-11260_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data3 <- read.table("SJ-2335-8553-nr4_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data4 <- read.table("UI-3077-A7551-20AT_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))

ACC6 <- as.data.table(rbind(data1, data2, data3, data4))
ACC6[ACC6$sample=="SJ-2335-11259-v1",]$sample <- "ACC6-P1"
ACC6[ACC6$sample=="UI-3077-11260",]$sample <- "ACC6-P2"
ACC6[ACC6$sample=="SJ-2335-8553-nr4",]$sample <- "ACC6-M1"
ACC6[ACC6$sample=="UI-3077-A7551-20AT",]$sample <- "ACC6-M2"

cnSpec(ACC6, y=custom)

#ACC7
data1 <- read.table("SJ-2335-11599_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("SJ-2335-2177-nr5_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC7 <- as.data.table(rbind(data1, data2))
ACC7[ACC7$sample=="SJ-2335-11599",]$sample <- "ACC7-P"
ACC7[ACC7$sample=="SJ-2335-2177-nr5",]$sample <- "ACC7-M"
cnSpec(ACC7, y=custom)

#ACC8
data1 <- read.table("SJ-2335-11691_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("UI-3077-11694_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data3 <- read.table("SJ-2335-11712_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data4 <- read.table("UI-3077-A9008-19_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data5 <- read.table("UI-3077-11912_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data6 <- read.table("WC-3639-11695_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC8 <- as.data.table(rbind(data1, data2, data3, data4, data5, data6))
ACC8[ACC8$sample=="SJ-2335-11691",]$sample <- "ACC8-P1"
ACC8[ACC8$sample=="UI-3077-11694",]$sample <- "ACC8-P2"
ACC8[ACC8$sample=="SJ-2335-11712",]$sample <- "ACC8-M1"
ACC8[ACC8$sample=="UI-3077-A9008-19",]$sample <- "ACC8-M2"
ACC8[ACC8$sample=="UI-3077-11912",]$sample <- "ACC8-M3"
ACC8[ACC8$sample=="WC-3639-11695",]$sample <- "ACC8-Thrombus"
cnSpec(ACC8, y=custom)


#ACC9
data1 <- read.table("UI-3077-11913_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data2 <- read.table("UI-3077-11913n2_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data3 <- read.table("UI-3077-11914_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
data4 <- read.table("UI-3077-11504_iCN.seg", header=F, col.names= c("sample", "chromosome", "start", "end", "dummy", "segmean"))
ACC9 <- as.data.table(rbind(data1, data2, data3, data4))
ACC9[ACC9$sample=="UI-3077-11913",]$sample <- "ACC9-P1"
ACC9[ACC9$sample=="UI-3077-11913n2",]$sample <- "ACC9-P2"
ACC9[ACC9$sample=="UI-3077-11914",]$sample <- "ACC9-R"
ACC9[ACC9$sample=="UI-3077-11504",]$sample <- "ACC9-M"
cnSpec(ACC9, y=custom)

cohort <- as.data.table(rbind(ACC1, ACC2, ACC3, ACC4, ACC5, ACC6, ACC7, ACC8, ACC9))
cnSpec(cohort, y = custom)
