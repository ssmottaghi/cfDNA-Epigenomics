Step 1: ATAC-seq Preprocessing
Script: 01_preprocess_bed.py

Input: Raw .bed files in Data/raw_bed/.

Action: Filters peaks for a score > 100 and removes duplicates based on genomic coordinates.

Output: Cleaned BED files with an n.bed suffix in Data/processed_bed/.



Step 2: Consensus Peak Generation
Script: 02_merge_peaks.sh

Requirement: This step requires the createIterativeOverlapPeakSet.R script from the [Corces Lab GitHub Repository](https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging/tree/main?tab=readme-ov-file).

Instructions:

Download createIterativeOverlapPeakSet.R from the link above.

Place it in the root directory of this project.

Run bash 02_merge_peaks.sh.

This will merge your individual preprocessed peaks into a non-overlapping master peak set.


Step 3: Signal Matrix Construction
Script: 03_build_count_matrix.R

Action: Standardizes all peaks to 501bp and overlaps individual sample peaks with the master consensus set.

Output: A numeric matrix (all_sam_mat.rds) where rows are genomic regions and columns are samples, containing the signal scores.



Step 4: Differential Accessibility Analysis

Script: 04_differential_analysis.R

Action: Performs DESeq2 normalization and runs "One-vs-All" contrasts for each tissue group. It identifies up to 2,000 unique, significant regions per tissue.

Output: A complete DESeq2 object and a filtered table of top markers (top_deg_median_expression.csv).


Step 5: Iterative Peak Selection & Heatmap
Script: 05_iterative_selection.R

Action: Implements a "rank and remove" strategy to select the top 200 unique genomic regions per tissue. This ensures no overlap between tissue-specific blocks.

Output: A publication-ready heatmap (Results/Tissue_Specific_Heatmap.pdf) and the final filtered peak matrix.


Step 6: cfDNA Depth Normalization
Script: 06_cfDNA_depth_normalization.R
Action: Calculates a normalized depth ratio for each ATAC-seq peak using a 5kb flanking window ($+/- 2500$ bp) in cfDNA samples. This identifies regions of differential nucleosome protection.
Input: CSV files in Data/cfDNA/ organized by offset.
Output: A normalized depth matrix (cfDNA_normalized_matrix.rds).


Step 7: Final Biomarker Cross-Validation
Script: 07_final_biomarker_selection.R

Action: Integrates ATAC-seq and cfDNA data. It selects the top 50 regions per tissue that show the most significant loss of nucleosome protection in cancer versus healthy cfDNA.
Output: The final biomarker matrix (Final_Biomarker_Matrix.rds) used for the paper's main conclusions.


Step 8: Reference Scaling & Deconvolution
Script: 08_cfDNA_deconvolution.R

Action: > 1. Creates a Reference Matrix (df_all_4) by Min-Max scaling the ATAC-seq signal of the final 650 regions.
2. Uses NNLS and SVR to deconvolve the patient cfDNA signal against this reference.
Output: Stacked bar plots in the Results/ folder showing the predicted tissue-of-origin for each cfDNA sample.


Step 9: Null Model Permutation Test
Script: 09_null_model_validation.R

Action: Runs 100 deconvolution iterations using random genomic regions. It compares the RCH (Relative Cancer/Healthy) values of the main model against this random distribution.
Significance: Proves that the identified biomarkers provide significantly higher signal-to-noise than random genomic background (represented by the "cfDECOR" star vs. the boxplots).
