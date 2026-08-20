#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=realSFS_het
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 24:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# per-sample heterozygosity from folded per-sample sfs

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
OUT="$PROJ/het"

cd "$OUT"

while read -r BAM; do
    SAMPLE=$(basename "$BAM" .bam)

    echo "[$(date)] Processing $SAMPLE"
    
    # angsd -i "$BAM" -ref "$REF" -anc "$REF" -P 16 \
    #   -GL 1 -doSaf 1 -doCounts 1 -minMapQ 20 -minQ 20 -setMinDepth 3 \
    #   -out "$OUT/$SAMPLE"

    realSFS -fold 1 "${SAMPLE}.saf.idx" > "${SAMPLE}_est.ml"

done < "$PROJ/bams_list.txt"
