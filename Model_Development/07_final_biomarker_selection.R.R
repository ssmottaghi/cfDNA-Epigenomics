# --- 1. Load Data ---
if(!exists("all_depth_df")) {
  all_depth_df <- readRDS("Data/Processed/cfDNA_normalized_matrix.rds")
  sign_matrix <- readRDS("Data/Processed/sign_matrix_selected.rds")
}

# 2. Compare Healthy (samples 1-3) vs Cancer (samples 4-8)
Normalied_depth_healthy <- c()
Normalied_depth_cancer <- c()
lenss <- c(0, cumsum(len_reg)) # len_reg was defined in Script 05

# Loop through cancer types (Prostate, Pancreas, Lung, Kidney, Colon)
for (i in 4:8){
  # Find which column in our matrix matches the cancer sample name
  j <- which(colnames(all_depth_df) == cfDNA_Samples[i])
  
  # Calculate average healthy signal for this specific tissue's 200 regions
  g1 <- rowMeans(all_depth_df[(lenss[j]+1): lenss[j+1], 1:3])
  # Get the signal for the cancer sample itself
  g2 <- all_depth_df[(lenss[j]+1): lenss[j+1], i]
  
  Normalied_depth_healthy <- c(Normalied_depth_healthy, g1)
  Normalied_depth_cancer <- c(Normalied_depth_cancer, g2)
}

# 3. Identify Top 50 regions based on cfDNA difference (Healthy - Cancer)
norm_depth_MN <- Normalied_depth_healthy - Normalied_depth_cancer
selected_region_names <- c()

for (i in 0:4) {
  # Take the block of 200 regions and narrow it down to the best 50
  highvec <- norm_depth_MN[(i*200 + 1): (i*200 + 200)]
  highvec_sorted <- sort(highvec, decreasing = TRUE)[1:50]
  selected_region_names <- c(selected_region_names, names(highvec_sorted))
}

# 4. Filter the original ATAC-seq signal matrix for these high-confidence regions
# sign_matrix_ordered should be created by subsetting your original normalized counts
final_selection_mat <- sign_matrix[selected_region_names, ]

# 5. Final Row-Normalization (Compositional Data)
# This scales the signal to proportions (0-1), making it ready for ML models
dffd_2 <- apply(final_selection_mat, 2, function(x) as.numeric(as.character(x)))
dffd_2 <- dffd_2 / (rowSums(dffd_2) + 1e-6) 

# --- 6. Final Export ---
if(!dir.exists("Data/Processed")) dir.create("Data/Processed")
saveRDS(dffd_2, "Data/Processed/Final_Biomarker_Matrix.rds")
write.csv(dffd_2, "Results/Final_Biomarkers_Proportions.csv")

print("Final biomarker selection complete. Ready for downstream modeling.")
