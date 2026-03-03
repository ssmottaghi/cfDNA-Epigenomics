import pysam
import numpy as np
import pandas as pd
import sys
import os
from collections import Counter

def run_analysis(bam_path, fasta_path, regions_file, depth_tag, out_prefix):
    samfile = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)
    
    # Storage for "One row per region" features
    epigenetic_data = []
    # Storage for "One row per base" WPS data
    wps_rows = []

    with open(regions_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3: continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            region_id = f"{chrom}_{start}_{end}"
            
            # Data for this specific region
            lengths, motifs, region_len = [], [], end - start
            spanning = np.zeros(region_len, dtype=int)
            endpoints = np.zeros(region_len, dtype=int)

            for read in samfile.fetch(chrom, max(0, start - 200), end + 200):
                if read.is_proper_pair and read.is_read1 and not read.is_duplicate:
                    f_start, tlen = read.reference_start, abs(read.template_length)
                    f_end = f_start + tlen
                    lengths.append(tlen)
                    
                    # End-Motif extraction
                    try: motifs.append(fasta.fetch(chrom, f_start, f_start + 4).upper())
                    except: pass

                    # WPS Calculation (Window=120)
                    for i in range(region_len):
                        pos = start + i
                        if f_start <= (pos - 60) and f_end >= (pos + 60): spanning[i] += 1
                        if ((pos - 60) <= f_start <= (pos + 60)) or ((pos - 60) <= f_end <= (pos + 60)):
                            endpoints[i] += 1

            # --- PROCESS WPS (Base-level) ---
            region_wps = spanning - endpoints
            for i, val in enumerate(region_wps):
                wps_rows.append({"chrom": chrom, "pos": start + i, "wps": val})

            # --- PROCESS EPIGENETICS (Region-level) ---
            if lengths:
                sfr = sum(1 for l in lengths if l < 100) / sum(1 for l in lengths if 150 <= l <= 220) if any(150 <= l <= 220 for l in lengths) else 0
                epigenetic_data.append({
                    "region": region_id,
                    "fragment_count": len(lengths),
                    "mean_len": np.mean(lengths),
                    "short_frag_ratio": sfr,
                    "occupancy": len([l for l in lengths if 140 <= l <= 200]) / region_len,
                    "top_motif": Counter(motifs).most_common(1)[0][0] if motifs else "N"
                })

    samfile.close()

    # Save Epigenetic Features CSV
    pd.DataFrame(epigenetic_data).to_csv(f"{out_prefix}_epigenetic_features.csv", index=False)
    # Save WPS CSV
    pd.DataFrame(wps_rows).to_csv(f"{out_prefix}_wps_base_level.csv", index=False)

if __name__ == "__main__":
    run_analysis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
