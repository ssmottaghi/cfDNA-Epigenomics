library(tidyverse)
library(nnls)
library(e1071)
library(RColorBrewer)
library(pheatmap)


# --- Preparation of the Reference Matrix (df_all_4) ---

# 1. Subset the original DESeq2 normalized counts (df) 
# using the final 650 regions selected in Step 7
# df contains the normalized counts for all samples and all regions
df_selected <- df[rownames(sign_matrix_based_on_cfDNA_ordered), ]

# 2. Convert to numeric matrix
# We ensure all columns are treated as numeric for the scaling math
df_all_2 <- apply(df_selected, 2, function(x) as.numeric(as.character(x)))

# 3. Min-Max Scaling and Inversion
# This step scales the ATAC-seq signal to a 0-1 range.
# The '1 - ...' inversion is used so that the 'highest' signal 
# in ATAC-seq (accessibility) matches the 'deepest' signal in cfDNA (protection loss).
df_all_4 <- 1 - ((df_all_2 - min(df_all_2)) / (max(df_all_2) - min(df_all_2)))

# Ensure the column names match your 13 Tissues/Cell-types
colnames(df_all_4) <- colnames(df_selected)
rownames(df_all_4) <- 1:nrow(df_all_4)
                  
# --- 2. Prepare Reference Matrix (A) and Mixture Matrix (b) ---
# df_all_4 represents your "Reference" (the ATAC-seq signal per tissue)


A <- apply(df_all_4, 2, function(x) as.numeric(as.character(x)))

# cfDNA signals to be deconvolved
all_depth_df_fin <- all_depth_df[rownames(sign_matrix_based_on_cfDNA_ordered),]
colnames(all_depth_df_fin) <- c("Healthy 1", "Healthy 2", "Healthy 3", "Prostate cancer", 
                               "Pancreas cancer", "Lung cancer", "Kidney cancer", "Colon cancer")

# Color Palette for the 13 tissues
my_colors_2 <- c("#2f3537", "#bfc1bb", "#00ffff", "#a40101", "#ef3636", "#ad7fa8", 
                 "#3566a4", "#4f9a07", "#8ae235", "#e9b96f", "#f57a01", "#c4a001","#fce950")

# --- 3. NNLS Deconvolution ---
message("Running NNLS Deconvolution...")
props_nnls <- data.frame(row.names = colnames(A))

for (col in 1:ncol(all_depth_df_fin)){
  b_raw <- all_depth_df_fin[, col]
  # Min-max scaling for the mixture vector
  b <- (b_raw - min(b_raw)) / (max(b_raw) - min(b_raw))
  
  result <- nnls(A, b)
  x <- result$x
  x <- x / sum(x) # Normalize to 1 (100%)
  props_nnls <- cbind(props_nnls, x)
}
colnames(props_nnls) <- colnames(all_depth_df_fin)

# --- 4. SVR Deconvolution ---
message("Running SVR Deconvolution...")
props_svr <- data.frame(row.names = colnames(A))

for (col in 1:ncol(all_depth_df_fin)){
  b <- all_depth_df_fin[, col]
  data_df <- data.frame(cbind(A, b))
  
  # Linear SVR model
  svr_model <- svm(b ~ ., data = data_df, type = "nu-regression", kernel = "linear")
  # Extract coefficients (slopes) for each tissue
  Prop <- coef(svr_model)[2:(ncol(A)+1)]
  Prop[Prop < 0] <- 0 # SVR can return negative values; we truncate to 0
  Prop <- Prop / sum(Prop)
  props_svr <- cbind(props_svr, Prop)
}
colnames(props_svr) <- colnames(all_depth_df_fin)

# --- 5. Plotting Function ---
# Creating a function to avoid repeating code for both plots
generate_stacked_bar <- function(props_mat, method_name, filename) {
  df_long <- props_mat %>%
    as.data.frame() %>%
    rownames_to_column("Tissue") %>%
    pivot_longer(-Tissue, names_to = "Sample", values_to = "Proportion")
  
  df_long$Sample <- factor(df_long$Sample, levels = colnames(props_mat))
  df_long$Tissue <- factor(df_long$Tissue, levels = rev(rownames(props_mat)))
  
  plt <- ggplot(df_long, aes(x = Sample, y = Proportion, fill = Tissue)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = rev(my_colors_2)) +
    labs(title = paste("Deconvolution by", method_name), 
         y = "Predicted Proportion", x = "cfDNA Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"))
  
  ggsave(filename, plt, width = 8, height = 6, dpi = 300)
}

# Generate final plots
generate_stacked_bar(props_nnls, "NNLS", "Results/Deconvolution_NNLS.tiff")
generate_stacked_bar(props_svr, "SVR", "Results/Deconvolution_SVR.tiff")
