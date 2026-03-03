library(ggfortify)
library(GenomicRanges)
library(readr)

# Helper function to read BED files into GRanges
readSummits2 <- function(file){
  # Columns: chr, start, end, name, score
  df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end", "name", "score"))))
  df <- df[,c(1,2,3,5)] # Exclude 'name' column to save memory
  return(GenomicRanges::makeGRangesFromDataFrame(df=df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE))
}

# 1. Load the Master Merged Peak set (Created by Step 02 Bash script)
# The master file is usually named 'All_Samples.fwp.filter.non_overlapping.bed' by the Corces script
merged_peak_path <- "Data/Results/All_Samples.fwp.filter.non_overlapping.bed"
merged_peak <- readSummits2(merged_peak_path)
merged_peak_d <- data.frame(merged_peak)

# 2. Load Metadata
file_info <- read.table('Data/metadata.txt', sep='\t', header=TRUE)
samples <- file_info$Sample

# 3. Initialize the Master Matrix
# Creating unique region identifiers (chr:start-end)
all_sam_mat <- data.frame(region=paste0(merged_peak_d$seqnames, ":", merged_peak_d$start, "-", merged_peak_d$end))
rownames(all_sam_mat) <- all_sam_mat$region

# 4. Loop through samples and fill the matrix using overlaps
# This maps the scores from individual 'n.bed' files to the master regions
for (i in samples) {
  # Path to the 'n.bed' files created in Step 01
  address <- paste0("Data/processed_bed/", i, "n.bed")
  
  if(file.exists(address)) {
    sample_peak <- readSummits2(address)
    # Standardize peak width to 501bp (250bp on each side of the summit)
    sample_peak <- resize(sample_peak, width = 501, fix = "center")
    
    # Identify which sample peaks overlap with which master peaks
    overlap <- data.frame(findOverlaps(sample_peak, merged_peak))
    
    new_col <- rep(0, length(merged_peak))
    sample_peak_d <- data.frame(sample_peak)
    
    # Assign the score from the sample to the corresponding master region
    new_col[overlap[,2]] = sample_peak_d$score[overlap[,1]]
    all_sam_mat[,i] <- new_col
  } else {
    warning(paste("File missing for sample:", i))
  }
}

# Remove the temporary 'region' column as they are now rownames
all_sam_mat$region <- NULL

# 5. Save the matrix
# We use .rds because it is compressed and preserves the R data structure perfectly
if(!dir.exists("Data/Processed")) dir.create("Data/Processed")
saveRDS(all_sam_mat, "Data/Processed/all_sam_mat.rds")

print("Signal matrix successfully created and saved to Data/Processed/all_sam_mat.rds")
