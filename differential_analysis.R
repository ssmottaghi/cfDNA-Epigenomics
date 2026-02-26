# ==============================================================================
# TISSUE-SPECIFIC FEATURE EXTRACTION (One-vs-All)
# ==============================================================================

# 1. Calculate median expression per group
cm_tissues <- t(aggregate(t(cm), list(file_info$Group), median))
colnames(cm_tissues) <- cm_tissues[1, ] # Fix header after aggregate
cm_tissues <- cm_tissues[-1, ] # Remove the group name row

# 2. Iterate through groups to find top differentially expressed genes
column_names = unique(file_info$Group) 
df = data.frame(matrix(nrow = 0, ncol = length(column_names))) 
colnames(df) = column_names
n = length(column_names)

for (i in 1:n) {
  # Create a contrast vector: 1 for current group, -1/(n-1) for all others
  vec <- rep((-1/(n-1)), n)
  vec[i] = 1
  
  # Run contrast
  ress = results(dds, lfcThreshold=1, contrast=vec)
  
  # Filter by LFC > 2 and padj < 0.05
  ress = ress[which(ress$log2FoldChange > 2 & ress$padj < 0.05), ]
  
  # Keep top 2000 genes by LFC
  ress = ress[rank(-abs(ress$log2FoldChange)) < 2001, ]
  
  # Add to summary table, avoiding duplicates
  if(nrow(ress) > 0) {
    duprows <- rownames(ress) %in% rownames(df)
    df <- rbind(df, cm_tissues[rownames(ress), ][!duprows, ])
  }
}
