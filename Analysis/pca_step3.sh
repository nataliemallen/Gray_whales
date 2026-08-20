#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=pca_ldprune
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# ld-pruned pca: ngsld + prune_graph, re-run pcangsd (Fig S2)

module --force purge
module load biocontainers
module load angsd
module load ngsld

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/PCA/LD_pruned"
mkdir -p "$OUT"; cd "$OUT"

# Number of individuals in the cohort
N_IND=$(wc -l < "$PROJ/bams_list.txt")

# Site positions from the genome-wide Beagle
zcat "$PROJ/PCA/whale.beagle.gz" | cut -f1 | tail -n +2 | sed 's/_/\t/' > whales.sites.txt
N_SITES=$(wc -l < whales.sites.txt)
# 
ngsLD --geno "$PROJ/PCA/whale.beagle.gz" --probs \
  --pos whales.sites.txt --max_kb_dist 10 --min_maf 0.05 --extend_out \
  --N_thresh 0.3 --call_thresh 0.9 --n_threads 60 --verbose 1 \
  --n_ind "$N_IND" --n_sites "$N_SITES" | sort -k1,1Vr -k2,2V > whales.allSNPs.ld

# Split contigs into large/small for prune_graph (memory)
awk '{print $1, $2}' "${REF}.fai" > contig_lengths.txt
awk '{if ($2 > 85000000) print $1 > "large_contigs.txt"; else print $1 > "small_contigs.txt"}' contig_lengths.txt

PRUNE=/home/allen715/prune_graph/target/release/prune_graph
while read -r contig; do
  grep -w "$contig" whales.allSNPs.ld > "${contig}.ld"
  $PRUNE --in "${contig}.ld" \
  --weight "column_6" \
  --filter "column_3 <= 10000 && column_6 >= 0.5" \
  --n-threads 40 \
  --out "${contig}_unlinked.ld"
done < large_contigs.txt
cat *_unlinked.ld > large_contigs_unlinked.ld

grep -Ff small_contigs.txt whales.allSNPs.ld > small_contigs.ld
$PRUNE --in small_contigs.ld \
  --weight "column_6" \
  --filter "column_3 <= 10000 && column_6 >= 0.5" \
  --n-threads 40 \
  --out small_contigs_unlinked.ld

cat large_contigs_unlinked.ld small_contigs_unlinked.ld > whales_unlinked.ld
sed 's/:/_/g' whales_unlinked.ld > whales_unlinked

# Subset the Beagle to unlinked sites
zcat "$PROJ/PCA/whale.beagle.gz" | head -n 1 > whales_unlinked.beagleheader
zcat "$PROJ/PCA/whale.beagle.gz" | grep -Fwf whales_unlinked > whales_unlinked.beagle
cat whales_unlinked.beagleheader whales_unlinked.beagle | gzip > whales_unlinked.beagle.gz

module load pcangsd
pcangsd -b whales_unlinked.beagle.gz -o whale_pruned -t 64

echo "Unlinked loci: $(wc -l < whales_unlinked)"
