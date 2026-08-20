#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=pcangsd
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 5-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# genome-wide covariance matrix from the beagle (Fig 4B)

module --force purge
module load biocontainers
module load angsd
module load pcangsd

PROJ="/scratch/gautschi/allen715/2026_whales"
cd "$PROJ/PCA"

pcangsd -b "$PROJ/PCA/whale.beagle.gz" -o whale -t 64
