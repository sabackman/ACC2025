library(ggplot2)
library(gplots)


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
activities <- read.table(file.path("spa_output", "Assignment_Solution", "Activities", "Assignment_Solution_Activities.txt"), header=T, row.names=1)

#Rename samples for publication ready figure
newnames <- rownames(activities)
for (i in 1:length(newnames)) {
  newnames[i] <- samplemapping[samplemapping$name==newnames[i],]$mapping
}
rownames(activities) <- newnames
normalized <- activities/rowSums(activities)
normalized <- normalized[,colSums(normalized)!=0] #Drop signatures without contribution

#All signatures with any contribution
#heatmap.2(as.matrix(normalized), trace="none", Rowv=FALSE, dendrogram="col", 
#          scale="none", col="bluered",
#          density.info="none", key.title=NA, margins=c(5,8))

#Only those with at least 10% contribution in at least one sample
heatmap.2(as.matrix(normalized[,colSums(normalized>=0.1)>0]), trace="none", Rowv=FALSE, dendrogram="col", 
          scale="none", col="bluered",
          density.info="none", key.title=NA, margins=c(5,8), key.xlab="Contribution (fraction)")
