#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=LD_prune_west
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

## prune LD

module purge
module load biocontainers
module load angsd
module load ngsld

# Set your input BAM file list
BAM_LIST="/scratch/negishi/allen715/Gray_whales/final_bams/merged/final_merged_use/new_merged_bams.list"
# Set the reference genome FASTA file path
REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
# Set output directory
OUT="/scratch/negishi/allen715/Gray_whales/PCA/LD_pruned/all/whales"

angsd -bam "$BAM_LIST" \
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
      -out "$OUT" \

zcat whales.beagle.gz | cut -f 1 | tail -n +2 | sed 's/_/\t/' > whales.pos.txt

#calculate ld for whales
east_IND=$(wc -l < east_bamlist.txt)
less whales.beagle.gz --pos | cut -f1 | sed 's/1_/1\t/g' | sed '1d' > whales.sites.txt
east_SITES=$(less whales.sites.txt | wc -l)

#if problem with extra line in sites file, do this, then redo site number calculation
head -n -1 whales.sites.txt > temp.txt ; mv temp.txt whales.sites.txt

ngsLD --geno whales.beagle.gz --probs \
    --pos whales.sites.txt \
    --max_kb_dist 10 --min_maf 0.05 --extend_out \
        --N_thresh 0.3 --call_thresh 0.9 \
        --n_threads 60 --verbose 1 \
        --n_ind 74 --n_sites 2465830 | \
             sort -k 1,1Vr -k 2,2V > whales.allSNPs.ld
            
# reference
REF="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"
# create contig_lengths.txt from .fai file
awk '{print $1, $2}' "${REF}.fai" > contig_lengths.txt
# extract contig names from whale.pos.txt
cut -f 1 whales.sites.txt | sort -u > contig_list.txt

# split contigs into "large" and "small"
awk '{if ($2 > 85000000) print $1 > "large_contigs.txt"; else print $1 > "small_contigs.txt"}' contig_lengths.txt

# process large contigs 
while read -r contig; do
    grep -w "$contig" whales.allSNPs.ld > "${contig}.ld"
    /home/allen715/prune_graph/target/release/prune_graph --in "${contig}.ld" \
        --weight-field "column_7" \
        --weight-filter "column_3 <= 10000 && column_7 >= 0.5" \
        --n-threads 40 \
        --out "${contig}_unlinked.ld"
done < large_contigs.txt

# combine large contig results
cat *_unlinked.ld > large_contigs_unlinked.ld

# process small contigs 
grep -wf small_contigs.txt whales.allSNPs.ld > small_contigs.ld
/home/allen715/prune_graph/target/release/prune_graph --in small_contigs.ld \
    --weight-field "column_7" \
    --weight-filter "column_3 <= 10000 && column_7 >= 0.5" \
    --n-threads 40 \
    --out small_contigs_unlinked.ld

# combine all results 
cat large_contigs_unlinked.ld small_contigs_unlinked.ld > whales_unlinked.ld

# make a copy 
cp whales_unlinked.ld whales_unlinked_copy.ld

# make output compatible with beagle
sed 's/:/_/g' whales_unlinked.ld > whales_unlinked

zcat whales.beagle.gz | head -n 1 > whales_unlinked.beagleheader  
zcat whales.beagle.gz | grep -Fwf whales_unlinked > whales_unlinked.beagle  
cat whales_unlinked.beagleheader whales_unlinked.beagle | gzip > whales_unlinked.beagle.gz  

module load pcangsd

pcangsd -b whales_unlinked.beagle.gz -o whale_pca -t 64
cat whales_unlinked | wc -l  # Count the number of unlinked loci
zcat whales_unlinked.beagle.gz | tail -n +2 | wc -l  # Check the number of rows (SNPs) in the final `beagle` file

