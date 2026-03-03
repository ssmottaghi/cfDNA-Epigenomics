# 1. Compare Healthy (samples 1-3) vs Cancer (samples 4-8)
Normalied_depth_healthy <- c()
Normalied_depth_cancer <- c()
lenss <- c(0, cumsum(len_reg))

# Loop through cancer types (e.g., Prostate, Pancreas, Lung, Kidney, Colon)
for (i in 4:8){
  j <- which(colnames(sign_matrix) == cfDNA_Samples[i])
  # Average healthy signal for these specific regions
  g1 <- rowMeans(all_depth_df[(lenss[j]+1): lenss[j+1], 1:3])
  # Specific cancer sample signal
  g2 <- all_depth_df[(lenss[j]+1): lenss[j+1], i]
  
  Normalied_depth_healthy <- c(Normalied_depth_healthy, g1)
  Normalied_depth_cancer <- c(Normalied_depth_cancer, g2)
}

# 2. Identify the Top 50 regions with the largest Healthy-Cancer difference
norm_depth_MN <- Normalied_depth_healthy - Normalied_depth_cancer
selected_regions <- c()

for (i in 0:4) {
  # Take the 200 regions from selection step 2 and pick the top 50 based on cfDNA
  highvec <- norm_depth_MN[(i*200 + 1): (i*200 + 200)]
  highvec <- sort(highvec, decreasing = TRUE)[1:50]
  selected_regions <- c(selected_regions, highvec)
}

# 3. Final Row-Normalization (Compositional Data)
# This converts absolute counts into relative proportions (0 to 1)
dffd_2 <- apply(sign_matrix_based_on_cfDNA_ordered, 2, function(x) as.numeric(as.character(x)))
dffd_2 <- dffd_2 / rowSums(dffd_2 + 1e-6) # Add epsilon to avoid division by zero

# Save the final biomarker matrix
saveRDS(dffd_2, "Data/Processed/Final_Biomarker_Matrix.rds")
