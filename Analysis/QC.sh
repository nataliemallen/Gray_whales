#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=aggregate_qc
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# cohort coverage summaries (breadth, depth) by population

module --force purge
module load biocontainers
module load samtools

PROJ="/scratch/gautschi/allen715/2026_whales"
QC="$PROJ/QC"
META="$PROJ/lists/sample_pop.tsv"
EXCLUDE="038046_ER-17-0237"

cd "$QC"

# Per-sample table straight from metadata (breadth/depth precomputed, mapping rate=1 placeholder
# unless you re-derive it from samtools flagstat below). Excludes the flagged sample.
# String-safe read loop so leading zeros (e.g. 015770) are never reformatted.
echo "sample,pop,breadth,depth" > qc_summary.csv
while IFS=$'\t' read -r sample pop breadth depth; do
  [ "$sample" = "sample" ] && continue
  [ -z "$sample" ] && continue
  [ "$sample" = "$EXCLUDE" ] && continue
  printf '%s,%s,%s,%s\n' "$sample" "$pop" "$breadth" "$depth" >> qc_summary.csv
done < "$META"

# OPTIONAL: regenerate breadth/depth directly from BAMs (uncomment to recompute from scratch).
# > qc_summary.csv ; echo "sample,pop,breadth,depth,maprate" > qc_summary.csv
# while IFS=$'\t' read -r s pop b d; do
#   [ "$s" = "sample" ] && continue
#   [ "$s" = "$EXCLUDE" ] && continue
#   bam="$PROJ/bams/${s}.bam"
#   # mean depth & breadth over covered genome
#   read depth breadth <<< $(samtools depth -a "$bam" | \
#       awk '{sum+=$3; if($3>0)cov++; tot++} END{printf "%.4f %.4f", sum/tot, 100*cov/tot}')
#   maprate=$(samtools flagstat "$bam" | awk '/mapped \(/{gsub(/[()%]/,"",$5); print $5; exit}')
#   echo "${s},${pop},${breadth},${depth},${maprate}" >> qc_summary.csv
# done < "$META"

# Population means + ranges
python3 - qc_summary.csv qc_summary_bypop.txt <<'PY'
import csv, sys, statistics as st
rows=list(csv.DictReader(open(sys.argv[1])))
def block(name, rs, fh):
    b=[float(r['breadth']) for r in rs]; d=[float(r['depth']) for r in rs]
    fh.write(f"{name}: n={len(rs)}\n")
    fh.write(f"  breadth  mean={st.mean(b):.2f}%  range={min(b):.2f}-{max(b):.2f}\n")
    fh.write(f"  depth    mean={st.mean(d):.2f}X  range={min(d):.2f}-{max(d):.2f}\n")
with open(sys.argv[2],'w') as fh:
    block("ALL", rows, fh)
    block("east (ENP)", [r for r in rows if r['pop']=='east'], fh)
    block("west (WNP)", [r for r in rows if r['pop']=='west'], fh)
print(open(sys.argv[2]).read())
PY
