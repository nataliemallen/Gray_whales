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

# for apparent parentage analysis

awk '{
    # Print the part of the line up to the "-9"
    for (i = 1; i <= 6; i++) printf "%s ", $i
    
    # Format the remaining columns with a slash between every pair
    for (i = 7; i <= NF; i += 2) {
        printf "%s/%s", $i, $(i + 1)
        if (i + 2 <= NF) printf " "
    }
    print ""
}' plink_output_copy.ped > apparent.txt
