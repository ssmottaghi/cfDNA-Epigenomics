#!/bin/bash
#SBATCH --job-name=cfDNA_ReadCounts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Match the environment we established for GCC 11.3
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0

# Path to your installed pysam
export PYTHONPATH=$PYTHONPATH:/mnt/isilon/cag_inf_op/palizbanf/misc/python_libs

ORIGINAL_BAM="/mnt/isilon/cag_inf_op/palizbanf/misc/IC37.bam"
REGIONS=(
    "/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all_regions 0.txt" 
    "/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all regions +2500.txt" 
    "/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all_regions -2500.txt"
)
DEPTHS=("original" "10" "5" "1" "0.5" "0.1")

for DEPTH in "${DEPTHS[@]}"; do
    if [ "$DEPTH" == "original" ]; then
        BAM_NAME="$ORIGINAL_BAM"
    else
        BAM_NAME="depth_${DEPTH}x.bam"
    fi

    # Skip if downsampled file hasn't been created yet
    if [ ! -f "$BAM_NAME" ]; then
        echo "BAM $BAM_NAME not found, skipping..."
        continue
    fi

    for REG_FILE in "${REGIONS[@]}"; do
        # Format filename to be safe (remove spaces)
        BASE_NAME=$(basename "$REG_FILE" .txt)
        SAFE_TAG=${BASE_NAME// /_}
        OUT_CSV="counts_depth_${DEPTH}x_${SAFE_TAG}.csv"
        
        echo "Processing $BAM_NAME for $REG_FILE"
        python calculate_depth_and_counts.py "$BAM_NAME" "$REG_FILE" "$OUT_CSV"
    done
done

echo "All counting tasks complete."
