# ==============================================================================
# SELECTION STEP 2: ITERATIVE RANKING AND DEDUPLICATION
# ==============================================================================

# cm_tissues_2: A copy of our median expression matrix to "drain" as we pick rows
cm_tissues_2 <- data.frame(cm_tissues) 
# df: The results from Selection Step 1 (top significant genes)
rownames(cm_tissues_2) <- rownames(df)

sign_matrix <- data.frame()
len_reg <- c()

# Loop through each tissue (n)
for (i in 1:n){
  # 1. Sort by highest expression in current tissue
  # We use the median expression values to find the most representative peaks
  new_rows <- cm_tissues_2[order(-cm_tissues_2[, i]), ][1:200, ]
  
  # 2. Append these 200 peaks to our final plotting matrix
  sign_matrix <- rbind(sign_matrix, new_rows)
  len_reg <- c(len_reg, nrow(new_rows))
  
  # 3. CRUCIAL: Remove these selected rows from the pool 
  # This prevents the same peak from appearing in two different tissue blocks
  om_ro <- which(rownames(cm_tissues_2) %in% rownames(new_rows))
  cm_tissues_2 <- cm_tissues_2[-om_ro, ]
}

# ==============================================================================
# FINAL PLOTTING
# ==============================================================================

pheatmap(sign_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("yellow", "white", "blue"))(100), 
         angle_col = 90, 
         show_rownames = FALSE,
         main = "Iterative Tissue-Specific Peaks (Top 200 per Group)")
