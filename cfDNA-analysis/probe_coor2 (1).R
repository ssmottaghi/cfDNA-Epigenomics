# manifest
library(readxl)

hm450.manifest.hg19<- read.csv("~/Downloads/humanmethylation450_15017482_v1-2.csv", header = T)
reference<-read.table("~/Desktop/projects/Panel_designer/mottaghi/reference.csv"
                      , sep="\033" , header = T)
ref_CpG<- reference$CpGs
# checking
length(ref_CpG)
length(which(ref_CpG%in%hm450.manifest.hg19$IlmnID))

out<- hm450.manifest.hg19[hm450.manifest.hg19$IlmnID%in% ref_CpG, ]
outdir= "~/Desktop/projects/Panel_designer/mottaghi/prob_coorhg19.txt"
#write.table(out, file = outdir, row.names = F , quote = F , sep = "\t")

#liftover to hg38
library(GenomicRanges)
d<- out[,c("CHR", "Coordinate_36", "Coordinate_36")]
rownames(d)<- out$IlmnID
colnames(d)<-c("seqnames", "start", "end")
d$seqnames<- paste0("chr", d$seqnames)
g<-makeGRangesFromDataFrame(d)

library(rtracklayer)
path<- "~/Desktop/hg19ToHg38.over.chain"
ch<- import.chain(path)

ldf<- liftOver(g,chain = ch)
outdir<-"~/Desktop/projects/Panel_designer/mottaghi/hg38.hm450.reference.txt"
write.table(as.data.frame(ldf), outdir , col.names = T , row.names = T
            , quote = F , sep = "\t")
