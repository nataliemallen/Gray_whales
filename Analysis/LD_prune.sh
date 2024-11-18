#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=LD_prune
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

### prune LD

module purge
module load biocontainers
module load angsd
module load ngsld

bam lists for east vs west
EAST="east_bamlist.txt"

EAST_OUT="/scratch/negishi/allen715/Gray_whales/PCA/LD_pruned/east"
REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"

# run angsd
angsd -bam "$EAST" \
      -doCounts 1 \
      -GL 1 \
      -doGlf 2 \
      -doMajorMinor 1 \
      -doMaf 1 \
      -minMaf 0.05 \
      -SNP_pval 1e-6 \
      -minMapQ 20 \
      -minQ 20 \
      -baq 1 \
      -ref "$REF" \
      -nThreads 60 \
      -out "$EAST_OUT" \

zcat east.beagle.gz | cut -f 1 | tail -n +2 | sed 's/_/\t/' > east.pos.txt

# calculate ld for EAST
east_IND=$(wc -l < east_bamlist.txt)
less east.beagle.gz --pos | cut -f1 | sed 's/1_/1\t/g' | sed '1d' > east.sites.txt
east_SITES=$(less east.sites.txt | wc -l)

# check for extra line in sites file, do this, then redo site number calculation
head -n -1 east.sites.txt > temp.txt ; mv temp.txt east.sites.txt

ngsLD --geno east.beagle.gz --probs \
    --pos east.sites.txt \
    --max_kb_dist 10 --min_maf 0.05 --extend_out \
        --N_thresh 0.3 --call_thresh 0.9 \
        --n_threads 60 --verbose 1 \
        --n_ind 33 --n_sites 2301096 | \
             sort -k 1,1Vr -k 2,2V > east.allSNPs.ld

# run prune_graph 
/home/allen715/prune_graph/target/release/prune_graph --in east.allSNPs.ld \
    --weight-field "column_7" \
    --weight-filter "column_3 <= 10000 && column_7 >= 0.5" \
    --out east_unlinked.ld

# make output compatible with beagle
sed 's/:/_/g' east_unlinked.ld > east_unlinked

zcat east.beagle.gz | head -n 1 > east_unlinked.beagleheader  # Extract the header
zcat east.beagle.gz | grep -Fwf east_unlinked > east_unlinked.beagle  # Keep only unlinked loci
cat east_unlinked.beagleheader east_unlinked.beagle | gzip > east_unlinked.beagle.gz  # Combine header and SNP data

pcangsd -b east_unlinked.beagle.gz -o east_pca -t ${N}
cat east_unlinked | wc -l  # Count the number of unlinked loci
zcat east_unlinked.beagle.gz | tail -n +2 | wc -l  # Check the number of rows (SNPs) in the final `beagle` file
