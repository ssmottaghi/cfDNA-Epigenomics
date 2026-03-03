#!/bin/bash
#SBATCH --job-name=cfDNA_3Reps
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --time=72:00:00

module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load Python/3.10.4-GCCcore-11.3.0

export PYTHONPATH=$PYTHONPATH:/mnt/isilon/cag_inf_op/palizbanf/misc/python_libs

HEALTHY="/mnt/isilon/cag_inf_op/palizbanf/misc/IH01.bam"
CANCER="/mnt/isilon/cag_inf_op/palizbanf/misc/IC37.bam"
REGIONS="/mnt/isilon/cag_inf_op/palizbanf/misc/mott/all_regions_0.txt"
FASTA="/mnt/isilon/cag_inf_op/palizbanf/misc/hs37d5.fa.gz"

TFs=("0.1" "0.05" "0.01" "0.005" "0.001")
REPS=("replicate_1" "replicate_2" "replicate_3")
SEEDS=("42" "77" "99")



for i in "${!REPS[@]}"; do
    REP=${REPS[$i]}
    SEED=${SEEDS[$i]}
    mkdir -p "$REP"

    for TF in "${TFs[@]}"; do
        OUT_BAM="$REP/sim_TF${TF}_${REP}.bam"
        
        # 1. Mix
        python simulate_mixing.py "$HEALTHY" "$CANCER" "$TF" 30 "$SEED" "$OUT_BAM"
        
        # 2. Verify Depth
        python calculate_depth_and_counts.py "$OUT_BAM" "$REGIONS" "$REP/depth_TF${TF}_${REP}.csv"
        
        # 3. Features
        python cfdna_engine.py "$OUT_BAM" "$FASTA" "$REGIONS" "30" "$REP/epi_TF${TF}_${REP}"
    done
done
