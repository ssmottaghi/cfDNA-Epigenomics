library(pheatmap)

# --- 1. Load Data from Step 4 ---
# We need the median expression matrix (cm_tissues) 
# and the group names (n) defined previously.
if(!exists("cm_tissues")) {
  load("Data/Processed/deseq_complete_results.RData")
}

cm_tissues_2 <- data.frame(cm_tissues) 
n <- ncol(cm_tissues_2)

sign_matrix <- data.frame()
len_reg <- c()

# --- 2. Iterative Ranking and Removal ---
for (i in 1:n){
  message(paste("Selecting top regions for:", colnames(cm_tissues_2)[i]))
  
  # Sort by highest median expression in the current tissue
  # We use [1:200, ] to take the top 200 peaks
  # na.last = NA handles cases where a tissue might have < 200 significant peaks
  ordered_rows <- cm_tissues_2[order(-cm_tissues_2[, i]), ]
  new_rows <- head(ordered_rows, 200)
  
  # Append to the final matrix
  sign_matrix <- rbind(sign_matrix, new_rows)
  len_reg <- c(len_reg, nrow(new_rows))
  
  # CRUCIAL: Remove these selected rows from the pool 
  # This prevents the same peak from appearing in two different tissue blocks
  om_ro <- which(rownames(cm_tissues_2) %in% rownames(new_rows))
  cm_tissues_2 <- cm_tissues_2[-om_ro, ]
}

# --- 3. Final Heatmap Generation ---
# Yellow-White-Blue scale is excellent for highlighting contrast in ATAC-seq signal
pheatmap(sign_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("yellow", "white", "blue"))(100), 
         angle_col = 90, 
         show_rownames = FALSE,
         main = "Iterative Tissue-Specific Peaks (Top 200 per Group)",
         filename = "Results/Tissue_Specific_Heatmap.pdf", # Optional: auto-save
         width = 8, height = 12)

# Save the selected peak names for the next step (cfDNA analysis)
saveRDS(sign_matrix, "Data/Processed/sign_matrix_selected.rds")
