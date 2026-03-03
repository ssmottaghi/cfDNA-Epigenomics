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
