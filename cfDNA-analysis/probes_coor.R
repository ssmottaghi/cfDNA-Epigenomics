# manifest
library(readxl)
load("~/Desktop/hm450.manifest.hg19.rda")
reference<-read.table("~/Desktop/projects/Panel_designer/mottaghi/reference.csv"
                       , sep="\033" , header = T)
ref_CpG<- reference$CpGs
# checking
length(ref_CpG)
length(which(ref_CpG%in%hm450.manifest.hg19$probeID))

out<- hm450.manifest.hg19[hm450.manifest.hg19$probeID%in% ref_CpG, ]
outdir= "~/Desktop/projects/Panel_designer/mottaghi/prob_coor.txt"
write.table(out, file = outdir, row.names = F , quote = F , sep = "\t")

#liftover to hg38
d<- out[,1:3]
colnames(d)<-c("seqnames", "start", "end")
g<-makeGRangesFromDataFrame(d)

path<- "~/Desktop/hg19ToHg38.over.chain"
ch<- import.chain(path)

ldf<- liftOver(g,chain = ch)
outdir<-"~/Desktop/projects/Panel_designer/mottaghi/hg38.hm450.manifest.txt"
write.table(as.data.frame(ldf), outdir , col.names = T , row.names = T
            , quote = F , sep = "\t")
