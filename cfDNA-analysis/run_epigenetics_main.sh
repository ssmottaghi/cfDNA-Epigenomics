#!/bin/bash
#SBATCH --job-name=cfDNA_Epi_Main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# 1. Load GCC 11.3 compatible modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0

# 2. Set Path for pysam
export PYTHONPATH=$PYTHONPATH:/mnt/isilon/cag_inf_op/palizbanf/misc/python_libs

ORIGINAL_BAM="/mnt/isilon/cag_inf_op/palizbanf/misc/IC37.bam"
BAM_10X="depth_10x.bam"
FASTA="/mnt/isilon/cag_inf_op/palizbanf/misc/hs37d5.fa.gz"

REGIONS=(
    "/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all_regions 0.txt" 
    "/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all regions +2500.txt" 
    "/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all_regions -2500.txt"
)

# List of BAMs to process
TARGET_BAMS=("$ORIGINAL_BAM" "$BAM_10X")

for BAM in "${TARGET_BAMS[@]}"; do
    # Create a label for the output (Original vs 10x)
    if [[ "$BAM" == *".bam" && "$BAM" != *"depth_"* ]]; then
        LABEL="Original"
    else
        LABEL="10x"
    fi

    for REG_FILE in "${REGIONS[@]}"; do
        # Extract filename and clean spaces
        BASE_NAME=$(basename "$REG_FILE" .txt)
        SAFE_TAG=${BASE_NAME// /_}
        
        OUT_PREFIX="epi_analysis_${LABEL}_${SAFE_TAG}"
        
        echo "================================================"
        echo "Running Full Epigenetics: $LABEL | $REG_FILE"
        echo "================================================"
        
        # Call the engine we built earlier
        python cfdna_engine.py "$BAM" "$FASTA" "$REG_FILE" "$LABEL" "$OUT_PREFIX"
    done
done
