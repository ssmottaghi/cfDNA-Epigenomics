import pysam
import numpy as np
import pandas as pd
import sys
import os

def get_region_metrics(bam_path, regions_file, out_csv):
    # Ensure BAM index exists
    if not os.path.exists(bam_path + ".bai"):
        print(f"Indexing {bam_path}...")
        pysam.index(bam_path)

    samfile = pysam.AlignmentFile(bam_path, "rb")
    results = []

    with open(regions_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3: continue
            
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            try:
                # 1. Simple Read Counting
                # .count() returns the number of alignments in the region
                read_count = samfile.count(chrom, start, end)

                # 2. Base-by-Base Depth Calculation
                counts = samfile.count_coverage(chrom, start, end)
                total_depth_per_base = np.sum(counts, axis=0)
                mean_depth = np.mean(total_depth_per_base) if len(total_depth_per_base) > 0 else 0
                
                results.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "read_count": read_count,
                    "avg_base_depth": round(mean_depth, 3),
                    "region_size": end - start
                })
            except ValueError:
                print(f"Region {chrom}:{start}-{end} not found in BAM header.")

    samfile.close()
    
    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv(out_csv, index=False)
    print(f"Saved results to {out_csv}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python calculate_depth_and_counts.py <bam> <regions> <output>")
    else:
        get_region_metrics(sys.argv[1], sys.argv[2], sys.argv[3])
