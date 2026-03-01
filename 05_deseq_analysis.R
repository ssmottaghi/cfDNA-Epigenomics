library(DESeq2)
library(ggfortify)
library(reshape2)

# 1. Load the matrix generated in Step 01
# (Assuming you saved it as an RDS file)
all_sam_mat <- readRDS("Data/Processed/all_sam_mat.rds")
file_info <- read.table('Data/metadata.txt', sep='\t', header=TRUE)

# 2. Pre-processing Metadata
file_info$Group <- gsub(" ", "", file_info$Group)
file_info$Group <- as.factor(file_info$Group)

# 3. Prepare the Matrix
# Add a pseudocount of 1 to avoid issues with zeros during modeling
all_sam_mat_3 <- as.matrix(all_sam_mat) + 1

# 4. Run DESeq2
dds = DESeqDataSetFromMatrix(countData = round(all_sam_mat_3), 
                              colData = file_info, 
                              design = ~ Group)

# We use fitType='local' as specified in our manuscript
dds = DESeq(dds, fitType='local')

# 5. Extract Normalized Counts
cm = data.frame(counts(dds, normalized=TRUE))

# Restore original zeros (remove the influence of the +1 pseudocount for visualization)
cm[all_sam_mat == 0] = 0


# Save the dds object and normalized counts for the plotting script
save(dds, cm, file = "Data/Processed/deseq_results.RData")


