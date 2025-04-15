#!/bin/bash
#SBATCH -A johnwayne
#SBATCH --job-name=whale_PLINK
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# Step 1: extract genotype data from .ped file
awk '{for(i=7;i<=NF;i+=2) printf "%s%s ", substr($i,1,1), substr($i+1,1,1); print ""}' plink_output.ped > temp_genotypes.txt

# Step 2: convert genotypes to COLONY format
awk '
{
    for(i=1;i<=NF;i++) {
        if($i == "00") printf "00 "
        else if($i == "11") printf "11 "
        else if($i == "12" || $i == "21") printf "12 "
        else if($i == "22") printf "22 "
        else printf "00 "  # For any unexpected values
    }
    print ""
}' temp_genotypes.txt > colony_genotypes.txt

# Step 3: add individual IDs and placeholder values for parents and sex
awk '{print NR, "0 0 0", $0}' colony_genotypes.txt > colony_input.txt

# Step 4: extract marker names from .map file
awk '{print $2}' plink_output.map > marker_names.txt

# Step 5: calculate allele frequencies
awk '
{
    for(i=1;i<=NF;i++) {
        if($i == "11") count[i"_1"] += 2
        else if($i == "12" || $i == "21") {count[i"_1"]++; count[i"_2"]++}
        else if($i == "22") count[i"_2"] += 2
        total[i] += 2
    }
}
END {
    for(i=1;i<=NF;i++) {
        freq1 = count[i"_1"] / total[i]
        freq2 = count[i"_2"] / total[i]
        if(freq1 == "") freq1 = 0
        if(freq2 == "") freq2 = 0
        print "1 2\n" freq1, freq2
    }
}' colony_genotypes.txt > allele_freqs.txt

# clean up temporary files
rm temp_genotypes.txt colony_genotypes.txt

echo "Conversion complete. Check colony_input.txt, marker_names.txt, and allele_freqs.txt"
