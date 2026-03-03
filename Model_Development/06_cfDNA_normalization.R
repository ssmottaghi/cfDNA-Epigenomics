library(tidyr)
library(readr)
library(dplyr)

# Helper function to read cfDNA depth files
readcfDNA <- function(file){
  if(!file.exists(file)) stop(paste("File not found:", file))
  # Using read_tsv with explicit col_types to ensure coordinates aren't mangled
  df <- suppressMessages(readr::read_tsv(file, col_names = FALSE))
  return(df)
}

# --- 1. Setup ---
# sign_matrix must be loaded from Step 5
if(!exists("sign_matrix")) {
  sign_matrix <- readRDS("Data/Processed/sign_matrix_selected.rds")
}

cfDNA_Samples <- c("Healthy", "IH02", "IH03", "Prostate", "Pancreas", "Lung", "Kidney", "Colon")
all_depth <- matrix(nrow = nrow(sign_matrix), ncol = length(cfDNA_Samples))
target_regions <- rownames(sign_matrix)

# --- 2. Iterative Normalization ---
for (j in seq_along(cfDNA_Samples)){
  sample_id <- cfDNA_Samples[j]
  message(paste("Processing cfDNA Sample:", sample_id))
  
  # --- Center Depth (0bp offset) ---
  # Expected format: "chr start end depth" OR "chr:start-end depth"
  path0 <- paste0("Data/cfDNA/0/all_in_", sample_id, ".csv")
  df0_raw <- readcfDNA(path0)
  colnames(df0_raw) <- c("chr_info", "depth")
  
  # Parse coordinates to standard chr1:start-end format
  df0_parsed <- df0_raw %>% 
    separate(chr_info, into = c("chr", "end"), sep = " (?=[^ ]+$)") %>%
    separate(chr, into = c("chr", "start"), sep = " (?=[^ ]+$)") %>%
    mutate(region_id = paste0("chr", chr, ":", start, "-", end))
  
  df0_vec <- setNames(df0_parsed$depth, df0_parsed$region_id)[target_regions]
  
  # --- Flank +2500bp ---
  path2500 <- paste0("Data/cfDNA/2500/all+2500_in_", sample_id, ".csv")
  df2500_raw <- readcfDNA(path2500)
  colnames(df2500_raw) <- c("chr", "start", "end", "depth")
  
  # Shift coordinates back to match the 'Center' labels
  df2500_vec <- df2500_raw %>%
    mutate(region_id = paste0("chr", chr, ":", start - 2500, "-", end - 2500)) %>%
    { setNames(.$depth, .$region_id) }[target_regions]
  
  # --- Flank -2500bp ---
  path_neg2500 <- paste0("Data/cfDNA/-2500/all-2500_in_", sample_id, ".csv")
  df_neg2500_raw <- readcfDNA(path_neg2500)
  colnames(df_neg2500_raw) <- c("chr", "start", "end", "depth")
  
  # Shift coordinates forward to match the 'Center' labels
  df_neg2500_vec <- df_neg2500_raw %>%
    mutate(region_id = paste0("chr", chr, ":", start + 2500, "-", end + 2500)) %>%
    { setNames(.$depth, .$region_id) }[target_regions]
  
  # --- 3. Normalization Formula ---
  # Ratio = (2 * Center) / (Flank_Left + Flank_Right)
  all_depth[,j] <- (2 * df0_vec) / (df2500_vec + df_neg2500_vec + 1e-6)
}

# --- 4. Final Formatting & Export ---
all_depth_df <- data.frame(all_depth)
colnames(all_depth_df) <- cfDNA_Samples
rownames(all_depth_df) <- target_regions

if(!dir.exists("Data/Processed")) dir.create("Data/Processed")
saveRDS(all_depth_df, "Data/Processed/cfDNA_normalized_matrix.rds")
