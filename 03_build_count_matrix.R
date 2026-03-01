library(ggfortify)
library(GenomicRanges)
library(readr)

# Helper function to read BED files into GRanges
readSummits2 <- function(file){
  df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end", "name", "score"))))
  df <- df[,c(1,2,3,5)] # Exclude 'name' column to save memory
  return(GenomicRanges::makeGRangesFromDataFrame(df=df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE))
}

# 1. Load the Master Merged Peak set
# Note: Path updated to be relative to the project root
merged_peak <- readSummits2("Data/Merged/All_Samples.fwp.filter.non_overlapping.bed")
merged_peak_d <- data.frame(merged_peak)

# 2. Load Metadata
file_info <- read.table('Data/metadata.txt', sep='\t', header=TRUE)
samples <- file_info$Sample
tissues <- file_info$Group

# 3. Initialize the Master Matrix
all_sam_mat <- data.frame(region=paste0(merged_peak_d$seqnames, ":", merged_peak_d$start, "-", merged_peak_d$end))
# Note: Use the regions as row names for easier indexing later
rownames(all_sam_mat) <- all_sam_mat$region

# 4. Loop through samples and fill the matrix using overlaps
for (i in samples) {
  address <- paste0("Data/Summits/", i, "n.bed")
  
  if(file.exists(address)) {
    sample_peak <- readSummits2(address)
    # Standardize peak width to 501bp (250bp each side + center)
    sample_peak <- resize(sample_peak, width = 501, fix = "center")
    
    # Calculate overlaps between sample peaks and master set
    overlap <- data.frame(findOverlaps(sample_peak, merged_peak))
    
    new_col <- rep(0, length(merged_peak))
    sample_peak_d <- data.frame(sample_peak)
    new_col[overlap[,2]] = sample_peak_d$score[overlap[,1]]
    all_sam_mat[,i] <- new_col
  } else {
    warning(paste("File missing for sample:", i))
  }
}

# Optional: Remove the temporary 'region' column if you set them as rownames
all_sam_mat$region <- NULL

# Save this matrix so the next script (DESeq2) can load it
saveRDS(all_sam_mat, "Data/all_sam_mat.rds")
