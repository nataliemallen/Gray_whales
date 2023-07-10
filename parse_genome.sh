#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

###parse reference genome file to check headers###

file="GCA_028021215.1_mEscRob2.pri_genomic.fa"

if [ ! -f "$file" ]; then
  echo "File not found: $file"
  exit 1
fi

grep "^>" "$file" | sed 's/^>//'
