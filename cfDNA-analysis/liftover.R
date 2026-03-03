# liftover
options(timeout = 1000)

library(rtracklayer)
library(liftOver)
library(GenomicRanges)

filenames<-list.files(path = "~/Desktop/projects/Panel_designer/mottaghi/mottaghi_files/")
setwd("~/Desktop/projects/Panel_designer/mottaghi/mottaghi_files/")
ldf <- lapply(filenames, read.table , header=FALSE)
ldf_modified <- lapply(ldf, function(df) {df$V1 <- paste0("chr", df$V1); df})
ldf_modified <- lapply(ldf_modified, function(df) {setNames(df, c("seqnames", "start", "end"))})

gr_list <- lapply(ldf_modified, function(df) makeGRangesFromDataFrame(df, keep.extra.columns = TRUE))


path<- "~/Desktop/hg19ToHg38.over.chain"
ch<- import.chain(path)

ldf_over<-lapply(gr_list , liftOver , chain=ch)

output_dir<-"./lifted/"
output_file <- paste0(output_dir, filenames)
for (i in seq_along(ldf_over)) {
  write.table(as.data.frame(gr_list[[i]]), file = output_file[i], sep = "\t", quote = FALSE, row.names = FALSE)
}
