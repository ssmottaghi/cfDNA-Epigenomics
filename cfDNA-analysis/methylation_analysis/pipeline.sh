#!/bin/bash

# Check if the number of input arguments is correct
if [ "$#" -ne 3 ]; then
   echo "executing!!!"
    exit 1
fi

acc=$1
r1=$2
r2=$3

mkdir $acc
BISMARK="/media/mehrmohammadi_hdd/mottaghi/Bismark-0.24.2/bismark"
REF="/media/second_hdd/Taklifi/genome_folder/hg38/fasta/default/"
DEDUP="/media/mehrmohammadi_hdd/mottaghi/Bismark-0.24.2/deduplicate_bismark"
METH="/media/mehrmohammadi_hdd/mottaghi/Bismark-0.24.2/bismark_methylation_extractor"
trimGalore="/media/mehrmohammadi_hdd/mottaghi/TrimGalore-0.6.10/trim_galore"

mkdir "${acc}/fastqc1"
# fastQC
fastqc --noextract -q -o "${acc}/fastqc1" $r1 
fastqc --noextract -q -o "${acc}/fastqc1" $r2

#echo  "fatsQC step completed"
#f_a="AGATTGGAAGAGTATATGTTT"
#r_a="AGATTGGAAGAGTGTTGTGTA"
#trimming
# default 
mkdir "${acc}/trimmed/"
#mkdir "${acc}/trimmed1/"
#$trimGalore --paired -o "${acc}/trimmed1/"  --basename $acc  $r1 $r2
#$trimGalore --paired -o "${acc}/trimmed1/"  --basename $acc -a $f_a -a2 $r_a  $r1 $r2
$trimGalore --paired --illumina  -o "${acc}/trimmed/"  --basename $acc  $r1 $r2

r1_t="${acc}/trimmed/${acc}_val_1.fq.gz"
r2_t="${acc}/trimmed/${acc}_val_2.fq.gz"


mkdir "${acc}/fastqc2"
fastqc  --noextract -q -o "${acc}/fastqc2"  $r1_t 
fastqc  --noextract -q -o "${acc}/fastqc2"  $r2_t 

# align to genome
align_out="${acc}/bam"
mkdir $align_out

$BISMARK  --genome $REF --multicore 20  -1 $r1_t -2 $r2_t  -o "${acc}/bam"
 
raw_bam="${acc}/bam/*_pe.bam"
# deduplicate
$DEDUP -p  --output_dir $align_out --bam $raw_bam

dedup_bam="${acc}/bam/*deduplicated.bam"
sorted_bam="${acc}/bam/${acc}.soted.bam"
soretd_bed="${acc}/bam/${acc}.soted.bed"

samtools sort-@ 20 $dedup_bam -o $sorted_bam
samtool index $sorted_bam 
bedtools bamtobed -i $sorted_bam > $sorted_bed
 
#extract methylations
mkdir "${acc}/meth"
$METH --gzip --bedGraph --ignore_r2 2 -o "${acc}/meth"  $dedup_bam 

bedGraph="${acc}/meth/${acc}_val_1_bismark_bt2_pe.deduplicated.bedGraph"

gunzip "${bedGraph}.gz"
sed -i '1d' $bedGraph

