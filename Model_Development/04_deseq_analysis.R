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

# fitType='local' is robust for ATAC-seq data where dispersion can vary
dds <- DESeq(dds, fitType = 'local')

# --- 3. Extract and Clean Normalized Counts ---
cm <- data.frame(counts(dds, normalized = TRUE))

# Restore original zeros to remove influence of pseudocount for visualization
cm[all_sam_mat == 0] <- 0

# --- 4. Calculate Median Expression per Group ---
# We use aggregate to find the median for each genomic region across tissue groups
cm_tissues <- t(aggregate(t(cm), list(file_info$Group), median))
colnames(cm_tissues) <- cm_tissues[1, ] 
cm_tissues <- cm_tissues[-1, ]          
# Convert to numeric matrix to prevent sorting/indexing errors
cm_tissues <- apply(cm_tissues, 2, as.numeric)
rownames(cm_tissues) <- rownames(cm)

# --- 5. Identify Top Differentially Expressed Regions ---
column_names <- levels(file_info$Group) 
n <- length(column_names)

# Initialize an empty data frame to store results
df_final_results <- data.frame(matrix(nrow = 0, ncol = n))
colnames(df_final_results) <- column_names

for (i in 1:n) {
  message(paste("Processing group:", column_names[i]))
  
  # One-vs-All contrast vector
  vec <- rep((-1/(n-1)), n)
  vec[i] <- 1
  
  # Run contrast
  ress <- results(dds, lfcThreshold = 1, contrast = vec)
  
  # Filter: LFC > 2 and padj < 0.05
  ress <- ress[which(ress$log2FoldChange > 2 & ress$padj < 0.05), ]
  
  if (nrow(ress) > 0) {
    # Keep top 2000 regions by Log2FoldChange
    ress <- ress[order(-ress$log2FoldChange), ]
    if (nrow(ress) > 2000) {
      ress <- ress[1:2000, ]
    }
    
    # Deduplicate: only add regions not already picked by a previous group
    new_regions <- setdiff(rownames(ress), rownames(df_final_results))
    
    if (length(new_regions) > 0) {
      df_final_results <- rbind(df_final_results, cm_tissues[new_regions, , drop = FALSE])
    }
  }
}

# --- 6. Save Final Outputs ---
if(!dir.exists("Data/Processed")) dir.create("Data/Processed")
save(dds, cm, df_final_results, file = "Data/Processed/deseq_complete_results.RData")
write.csv(df_final_results, "Data/Processed/top_deg_median_expression.csv")

print("Differential analysis complete. Results saved to Data/Processed/")
