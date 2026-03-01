import sys
import os
import subprocess

def run_cmd(cmd):
    print(f"Executing: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def mix_samples(healthy_bam, cancer_bam, tf, target_depth, seed, out_bam):
    # Current depth estimates (Based on your previous idxstats output)
    h_curr_cov = 30.0  
    c_curr_cov = 45.0

    # Calculate required fractions
    h_frac = ((1 - tf) * target_depth) / h_curr_cov
    c_frac = (tf * target_depth) / c_curr_cov

    # Samtools -s SEED.FRACTION
    h_arg = f"{seed}{str(round(h_frac, 4))[1:]}"
    c_arg = f"{seed + 10}{str(round(c_frac, 4))[1:]}"

    h_tmp = f"h_tmp_{tf}_{seed}.bam"
    c_tmp = f"c_tmp_{tf}_{seed}.bam"

    # Subsample and Merge
    run_cmd(f"samtools view -b -s {h_arg} -o {h_tmp} {healthy_bam}")
    run_cmd(f"samtools view -b -s {c_arg} -o {c_tmp} {cancer_bam}")
    run_cmd(f"samtools merge -f {out_bam} {h_tmp} {c_tmp}")
    run_cmd(f"samtools index {out_bam}")

    # Cleanup
    if os.path.exists(h_tmp): os.remove(h_tmp)
    if os.path.exists(c_tmp): os.remove(c_tmp)

if __name__ == "__main__":
    # Args: healthy_bam, cancer_bam, tf, target_depth, seed, out_bam
    mix_samples(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6])
