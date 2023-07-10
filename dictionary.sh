#!/bin/bash
#SBATCH --job-name=whale_dictionary
#SBATCH -A fnrpredator
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH --error=dict.err
#SBATCH --output=dict.out
#SBATCH --job-name=produce_align_SLURMM_jobs
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load bioinfo
module load bwa
module load picard-tools
module load samtools

REF=/scratch/bell/allen715/Gray_whales/Reference/original.fa

bwa index $REF
samtools faidx $REF
PicardCommandLine CreateSequenceDictionary reference=$REF output=original.fa.dict
