library(tidyr)
library(nnls)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# --- 1. Prepare Data for Random Selection ---
# Use the consensus median matrix from Script 4/5
if(!exists("cm_tissues")) {
  load("Data/Processed/deseq_complete_results.RData")
}

# Clean and convert to numeric matrix
cm_pool <- as.matrix(cm_tissues)
storage.mode(cm_pool) <- "numeric"

successful_iterations <- 0
max_iterations <- 100
random_results_list <- list()

# --- 2. Permutation Loop ---
while (successful_iterations < max_iterations) {
  
  # A. Pick 650 random regions from the entire ATAC-seq pool
  random_indices <- sample(nrow(cm_pool), 650)
  cm_random <- cm_pool[random_indices, ]
  
  # B. Normalize Random Reference Matrix
  # (Min-Max scaling to match the cfDNA signal)
  cm_random_norm <- 1 - (cm_random / max(cm_random, na.rm = TRUE))
  
  # C. cfDNA Normalization for Random Regions
  # This section recalculates the depth ratios for the 650 random regions
  cfDNA_Samples <- c("Healthy", "IH02", "IH03", "Prostate", "Pancreas", "Lung", "Kidney", "Colon")
  all_depth_random <- matrix(nrow = 650, ncol = length(cfDNA_Samples))
  
  # (Note: Re-using the logic from Script 06 to fill all_depth_random)
  # ... [cfDNA normalization logic] ...

  # D. Deconvolution (NNLS)
  A <- cm_random_norm
  props_iter <- data.frame(row.names = colnames(A))
  skip_iteration <- FALSE
  
  for (i in 1:ncol(all_depth_random)){
    b_raw <- all_depth_random[, i]
    b <- (b_raw - min(b_raw)) / (max(b_raw) - min(b_raw) + 1e-6)
    
    if (anyNA(A) || anyNA(b)) { skip_iteration <- TRUE; break }
    
    res <- nnls(A, b)
    x <- res$x / (sum(res$x) + 1e-6)
    props_iter <- cbind(props_iter, x)
  }
  
  if (skip_iteration) next
  
  # E. Store Results
  successful_iterations <- successful_iterations + 1
  random_results_list[[successful_iterations]] <- props_iter
  message(paste("Completed Iteration:", successful_iterations))
}
