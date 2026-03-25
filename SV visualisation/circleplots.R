library(circlize)
library(data.table)
source("plotfunction.R")

samples <- list.files("merged_sv_vcfs")

samplenames <- c("SJ-2335-7504", "UI-3077-7510", "SJ-2335-8745", "SJ-2335-8748", "SJ-2335-8749", 
                 "SJ-2335-10817-2", "SJ-2335-10820-2",
                 "SJ-2335-5589", "SJ-2335-6031", 
                 "SJ-2335-2269", "SJ-2335-4465", 
                 "SJ-2335-18381-nr3", "SJ-2335-11529",
                 "SJ-2335-11259", "UI-3077-11260", "SJ-2335-8553-nr4", "UI-3077-A7551-20AT", 
                 "SJ-2335-11599", "SJ-2335-2177-nr5", 
                 "SJ-2335-11691", "UI-3077-11694", "SJ-2335-11712", "UI-3077-A9008-19", "UI-3077-11912", "WC-3639-11695",
                 "UI-3077-11913", "UI-3077-11913n2", "UI-3077-11914", "UI-3077-11504")
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


for (s in samples) {
  patient <- strsplit(s, "\\.")[[1]][1]
  samplename <- gsub(".WC-3639-", "", gsub("UI-3077-", "", gsub("SJ-2335-", "", strsplit(s, "\\.")[[1]][2])))
  samplename_ <- gsub(".vcf", "", strsplit(s, "\\.")[[1]][2])
  message(samplename, " ", samplemapping[samplemapping$name==samplename_,]$mapping)
  png(file.path("plots", paste(paste(patient, samplename, sep="-"), ".png", sep="")), res=300, width=15, height=15, units="cm")
  circos.par(circle.margin=0.05)
  plot_genome(file.path("merged_sv_vcfs", s), samplemapping[samplemapping$name==samplename_,]$mapping)
  dev.off()
}
