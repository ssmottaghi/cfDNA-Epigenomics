library(GenomicRanges)
library(readr)

# Function to read summits into GRanges
readSummits2 <- function(file){
  df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end", "name", "score"))))
  df <- df[,c(1,2,3,5)] # Drop 'name' to save memory
  return(GenomicRanges::makeGRangesFromDataFrame(df=df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE))
}

# 1. Load the merged peak set (Result from the Corces script)
# Update this path to match your repo structure
merged_peak_path <- "Data/ATACdb/Results/All_Samples.fwp.filter.non_overlapping.bed"
merged_peak <- readSummits2(merged_peak_path)
merged_peak_d <- data.frame(merged_peak)

# 2. Load Metadata
file_info <- read.table('Data/ATACdb/metadata_new_3.txt', sep='\t', header=TRUE)
samples <- unique(file_info$Group)

# 3. Build the Signal Matrix
sign_mat <- data.frame(region=paste0(merged_peak_d$seqnames, ":", merged_peak_d$start, "-", merged_peak_d$end))

for (i in samples) {
  # Construct address based on sample name
  address <- paste0("Data/ATACdb/Results/", i, "/", i, ".fwp.filter.non_overlapping.bed")
  
  if(file.exists(address)) {
    sample_peak <- readSummits2(address)
    sample_peak <- resize(sample_peak, width = 2 * 250 + 1, fix = "center") # Standardize width
    
    overlap <- data.frame(findOverlaps(sample_peak, merged_peak))
    
    new_col <- rep(0, length(merged_peak))
    sample_peak_d <- data.frame(sample_peak)
    new_col[overlap[,2]] = sample_peak_d$score[overlap[,1]]
    sign_mat[,i] <- new_col
  } else {
    warning(paste("File not found for sample:", i))
  }
}
