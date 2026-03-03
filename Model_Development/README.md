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
