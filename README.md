Step 1: ATAC-seq Preprocessing
Script: 01_preprocess_bed.py

Input: Raw .bed files in Data/raw_bed/.

Action: Filters peaks for a score > 100 and removes duplicates based on genomic coordinates.

Output: Cleaned BED files with an n.bed suffix in Data/processed_bed/.




Dependencies
1. For peak merging:
https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging
Download createIterativeOverlapPeakSet.R from https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging and place it in the /tools folder
