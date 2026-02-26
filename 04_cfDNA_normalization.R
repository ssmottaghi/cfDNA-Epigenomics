library(tidyr)
library(readr)

# Function to read cfDNA depth files
readcfDNA <- function(file){
  # Ensure the file exists before reading
  if(!file.exists(file)) stop(paste("File not found:", file))
  df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr_info", "depth"))))
  return(df)
}

# 1. Setup Samples and Containers
cfDNA_Samples <- c("Healthy", "IH02", "IH03", "Prostate", "Pancreas", "Lung", "Kidney", "Colon")
# sign_matrix must be loaded from the previous script
all_depth <- matrix(nrow = nrow(sign_matrix), ncol = length(cfDNA_Samples))

# 2. Loop through each sample to calculate Normalized Depth
for (j in seq_along(cfDNA_Samples)){
  i <- cfDNA_Samples[j]
  
  # --- Center Depth (0) ---
  path0 <- paste0("Data/cfDNA/0/all_in_", i, ".csv")
  df2_raw <- readcfDNA(path0)
  # Parsing the coordinate string (e.g., "1 100 200") to match sign_matrix rownames
  df2_split <- df2_raw %>% 
    separate(chr_info, into = c("chr", "end"), sep = " (?=[^ ]+$)") %>%
    separate(chr, into = c("chr", "start"), sep = " (?=[^ ]+$)")
  
  names2 <- paste0("chr", df2_split$chr, ":", df2_split$start, "-", df2_split$end)
  df2 <- setNames(df2_raw$depth, names2)[rownames(sign_matrix)]
  
  # --- Flank +2500 ---
  path2500 <- paste0("Data/cfDNA/2500/all+2500_in_", i, ".csv")
  df2500_raw <- readcfDNA(path2500)
  # Shift coordinates back to center for matching
  names2500 <- paste0("chr", df2500_raw$X1, ":", df2500_raw$X2 - 2500, "-", df2500_raw$X3 - 2500)
  df2500 <- setNames(df2500_raw$depth, names2500)[rownames(sign_matrix)]
  
  # --- Flank -2500 ---
  path_neg2500 <- paste0("Data/cfDNA/-2500/all-2500_in_", i, ".csv")
  df_neg2500_raw <- readcfDNA(path_neg2500)
  # Shift coordinates back to center for matching
  names_neg2500 <- paste0("chr", df_neg2500_raw$X1, ":", df_neg2500_raw$X2 + 2500, "-", df_neg2500_raw$X3 + 2500)
  df_neg2500 <- setNames(df_neg2500_raw$depth, names_neg2500)[rownames(sign_matrix)]
  
  # 3. Calculate Normalized Depth: (2 * Center) / (Flank1 + Flank2)
  # Adding a small epsilon (1e-6) to denominator avoids division by zero
  all_depth[,j] <- (2 * df2) / (df_neg2500 + df2500 + 1e-6)
}

# 4. Final Formatting
all_depth_df <- data.frame(all_depth)
colnames(all_depth_df) <- cfDNA_Samples
rownames(all_depth_df) <- rownames(sign_matrix)

# Save for downstream plotting or machine learning
saveRDS(all_depth_df, "Data/Processed/cfDNA_normalized_matrix.rds")
