#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=outgroup_map
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# map fin whale (SRR23615109) to the gray whale reference as phylogeny outgroup

module --force purge
module load biocontainers
module load sra-tools
module load fastqc
module load bwa
module load samtools

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/reference/outgroup"
SRA="SRR23615109"
THREADS=64
mkdir -p "$OUT"; cd "$OUT"

# 1) download + extract paired FASTQs
# prefetch --max-size u "$SRA"
# fasterq-dump --split-files -e "$THREADS" -O . "$SRA"

# 2) read QC
# fastqc -t "$THREADS" "${SRA}_1.fastq" "${SRA}_2.fastq"

# 3) map, sort, mark duplicates
# bwa mem -t "$THREADS" -R "@RG\tID:${SRA}\tSM:finwhale\tPL:ILLUMINA" \
#   "$REF" "${SRA}_1.fastq" "${SRA}_2.fastq" \
#   | samtools sort -@ "$THREADS" -o ${SRA}.sorted.bam -

# name-sort
samtools sort -n -@ "$THREADS" \
    -o ${SRA}.namesort.bam \
    ${SRA}.sorted.bam

# add mate information required by markdup
samtools fixmate -m \
    ${SRA}.namesort.bam \
    ${SRA}.fixmate.bam

# coordinate sort
samtools sort -@ "$THREADS" \
    -o ${SRA}.sorted2.bam \
    ${SRA}.fixmate.bam

# mark duplicates
samtools markdup -@ "$THREADS" \
    ${SRA}.sorted2.bam \
    ${SRA}.markdup.bam

samtools index ${SRA}.markdup.bam

# remove intermediate BAMs to save space
rm ${SRA}.namesort.bam ${SRA}.fixmate.bam ${SRA}.sorted2.bam

# 4 remove unmapped/secondary/QCfail/duplicates/supplementary
samtools view -@ "$THREADS" \
    -b \
    -f 0x2 \
    -F 0x904 \
    ${SRA}.markdup.bam \
    > finwhale.bam

samtools index finwhale.bam
