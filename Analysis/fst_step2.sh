#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=fst_windows
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 06:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# bin per-site fst into 500 kb windows for jackknife/subsampling significance (Table S2)

PROJ="/scratch/gautschi/allen715/2026_whales"
cd "$PROJ/fst"

# realSFS fst print columns: chr  pos  A(alpha)  B(beta)
awk 'BEGIN{OFS="\t"} {
  w=int($2/500000)*500000; key=$1"\t"w
  A[key]+=$3; B[key]+=$4; N[key]++
} END{ for(k in A) print k, A[k], B[k], N[k] }' east_west_fst_persite.txt \
  | sort -k1,1V -k2,2n > windows_AB.txt

echo "500kb windows: $(wc -l < windows_AB.txt)"
echo "Global weighted FST (sumA/sumB): $(awk '{a+=$3;b+=$4} END{print a/b}' windows_AB.txt)"

