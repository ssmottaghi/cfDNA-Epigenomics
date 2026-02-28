#!/bin/bash

# 1. Ensure the directories match your project structure
# We point to the output of Script 01 (processed_bed) as the macs2dir
METADATA="Data/metadata.txt"
SUMMIT_DIR="Data/processed_bed/"
OUT_DIR="Data/ATACdb/Results/"
BLACKLIST="Data/blacklist19.bed"

# 2. Check if the external R script is present
if [ ! -f "createIterativeOverlapPeakSet.R" ]; then
    echo "Error: createIterativeOverlapPeakSet.R not found in the current directory."
    exit 1
fi

# 3. Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# 4. Run the iterative overlap script
# Note: --suffix "n.bed" matches the output from your Python script
Rscript createIterativeOverlapPeakSet.R \
  --metadata $METADATA \
  --macs2dir $SUMMIT_DIR \
  --outdir $OUT_DIR \
  --suffix "n.bed" \
  --blacklist $BLACKLIST \
  --genome "hg19" \
  --spm 5 \
  --rule "(n+1)/2" \
  --extend 250
