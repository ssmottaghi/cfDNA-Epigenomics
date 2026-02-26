Rscript createIterativeOverlapPeakSet.R \
  --metadata Data/ATACdb/metadata.txt \
  --macs2dir Data/ATACdb/Summits/ \
  --outdir Data/ATACdb/Results/ \
  --suffix "n.bed" \
  --blacklist Data/ATACdb/blacklist19.bed \
  --genome "hg19" \
  --spm 5 \
  --rule "(n+1)/2" \
  --extend 250
