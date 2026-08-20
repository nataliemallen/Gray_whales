#!/bin/bash
# qc step 2: cohort relatedness, inbreeding, duplicates, then aggregate to a master table
#SBATCH -A dewoody
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -c 128
#SBATCH -t 1-0:00:00
#SBATCH --job-name=qc2_cohort
#SBATCH -e logs/%x_%j.err
#SBATCH -o logs/%x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

set -uo pipefail
module --force purge
module load biocontainers
module load angsd
module load pcangsd
module load samtools
module load bcftools

module load anaconda
conda activate ngsrelate_env 

# module load samtools angsd
mkdir -p logs

# config
bam_list="/scratch/gautschi/allen715/whale/bams_list.txt"      # one BAM path per line
ref_fasta="/scratch/gautschi/allen715/2026_whales/reference/ref.fa"
outdir="/scratch/gautschi/allen715/2026_whales/QC"
pop_csv="/scratch/gautschi/allen715/2026_whales/whale_nobo_pop.csv"   # ind,sample_ID,pop (same order as bam_list!)
THREADS=${SLURM_CPUS_PER_TASK:-128}
MINMAPQ=30
MINQ=20
N_IND=$(wc -l < "$bam_list")     # <-- number of BAMs == ngsRelate -n. DO NOT hardcode.
echo "N individuals = $N_IND"

mkdir -p "$outdir/cohort"
cd "$outdir/cohort"

# genotype likelihoods
# run A: glf3 + mafs for ngsRelate
angsd -bam "$bam_list" -ref "$ref_fasta" -GL 2 \
      -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
      -doGlf 3 -minMapQ $MINMAPQ -minQ $MINQ -nThreads $THREADS -out cohort_glf3
zcat cohort_glf3.mafs.gz | cut -f6 | sed 1d > freq

# run B: beagle (glf2) for PCAngsd
angsd -bam "$bam_list" -ref "$ref_fasta" -GL 2 \
      -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
      -doGlf 2 -minMapQ $MINMAPQ -minQ $MINQ -nThreads $THREADS -out cohort_beagle

# (5) relatedness
ngsRelate -g cohort_glf3.glf.gz -n "$N_IND" -f freq -O relate.res -p $THREADS

# (6) inbreeding f + pca (pcangsd)
# (pcangsd flag names vary by version; --inbreedSamples writes per-sample F)
pcangsd -b cohort_beagle.beagle.gz -o pcangsd -t $THREADS --inbreedSamples || \
pcangsd --beagle cohort_beagle.beagle.gz --out pcangsd --threads $THREADS --inbreed-samples || \
echo "WARN: pcangsd inbreeding step failed - check flag syntax for your version"

# (7) cohort ibs matrix for duplicate/swap detection
angsd -bam "$bam_list" -ref "$ref_fasta" -GL 1 \
      -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
      -doIBS 1 -makeMatrix 1 -doCounts 1 \
      -minMapQ $MINMAPQ -minQ $MINQ -nThreads $THREADS -out cohort_ibs

# aggregate everything
python3 "/scratch/gautschi/allen715/2026_whales/aggregate_qc.py" \
  --persample_dir "$outdir/persample" \
  --relate relate.res \
  --pop "$pop_csv" \
  --inbreed pcangsd.inbreed.samples \
  --ibs cohort_ibs.ibsMat \
  --bamlist "$bam_list" \
  --out "$outdir/MASTER_QC_SUMMARY.tsv" \
  --flags "$outdir/FLAGGED_SAMPLES.txt"

echo ""
echo "================================================================"
echo "MASTER summary : $outdir/MASTER_QC_SUMMARY.tsv"
echo "Flagged report : $outdir/FLAGGED_SAMPLES.txt"
echo "================================================================"
cat "$outdir/FLAGGED_SAMPLES.txt"

