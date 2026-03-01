library(DESeq2)
library(ggfortify)
library(reshape2)

# --- 1. Data Loading & Pre-processing ---
all_sam_mat <- readRDS("Data/Processed/all_sam_mat.rds")
file_info <- read.table('Data/metadata.txt', sep='\t', header=TRUE)

# Clean up Group names and ensure it is a factor
file_info$Group <- as.factor(gsub(" ", "", file_info$Group))

# Add pseudocount of 1 to handle zeros during modeling
all_sam_mat_3 <- as.matrix(all_sam_mat) + 1

# --- 2. Run DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = round(all_sam_mat_3), 
                              colData = file_info, 
                              design = ~ Group)

# fitType='local' as specified in your manuscript
dds <- DESeq(dds, fitType = 'local')

# --- 3. Extract and Clean Normalized Counts ---
cm <- data.frame(counts(dds, normalized = TRUE))

# Restore original zeros to remove influence of pseudocount for visualization
cm[all_sam_mat == 0] <- 0

# --- 4. Calculate Median Expression per Group ---
# We transpose 'cm' to aggregate by Group, then transpose back
cm_tissues <- t(aggregate(t(cm), list(file_info$Group), median))
colnames(cm_tissues) <- cm_tissues[1, ] # Set group names as headers
cm_tissues <- cm_tissues[-1, ]          # Remove the temporary group name row
cm_tissues <- as.matrix(cm_tissues)     # Convert to matrix for easier indexing

# --- 5. Identify Top Differentially Expressed Genes (DEG) ---
column_names <- unique(file_info$Group) 
n <- length(column_names)

# Initialize an empty data frame to store results
df_final_results <- data.frame(matrix(nrow = 0, ncol = n))
colnames(df_final_results) <- column_names

for (i in 1:n) {
  # Create a contrast vector: 1 for current group vs average of others
  vec <- rep((-1/(n-1)), n)
  vec[i] <- 1
  
  # Run contrast with LFC threshold
  ress <- results(dds, lfcThreshold = 1, contrast = vec)
  
  # Filter by Significance: LFC > 2 and padj < 0.05
  ress <- ress[which(ress$log2FoldChange > 2 & ress$padj < 0.05), ]
  
  # Keep top 2000 genes by Log2FoldChange
  if (nrow(ress) > 0) {
    ress <- ress[order(-ress$log2FoldChange), ]
    if (nrow(ress) > 2000) {
      ress <- ress[1:2000, ]
    }
    
    # Identify genes not already in our final table to avoid duplicates
    new_genes <- setdiff(rownames(ress), rownames(df_final_results))
    
    # Add median expression data for these new genes to the final table
    if (length(new_genes) > 0) {
      df_final_results <- rbind(df_final_results, cm_tissues[new_genes, , drop = FALSE])
    }
  }
}

# --- 6. Save Final Outputs ---
save(dds, cm, df_final_results, file = "Data/Processed/deseq_complete_results.RData")
write.csv(df_final_results, "Data/Processed/top_deg_median_expression.csv")





