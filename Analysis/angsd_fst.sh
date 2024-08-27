#!/bin/bash
#SBATCH --job-name=fst
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --error=fst.err
#SBATCH --output=fst.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load biocontainers
module load angsd

#!/bin/bash

# List of site identifiers
sites=("east" "west")

# Path to reference genome
reference="/scratch/negishi/allen715/Gray_whales/reference/ref.fa"

# Log file to track progress
log_file="progress.log"

# Function to check if a step is complete
is_step_complete() {
    grep -q "$1" "$log_file"
}

# Function to log a completed step
log_step() {
    echo "$1" >> "$log_file"
}

# Loop through each pair of sites
for ((i=0; i<${#sites[@]}; i++)); do
    for ((j=i+1; j<${#sites[@]}; j++)); do
        site1=${sites[$i]}
        site2=${sites[$j]}

        # Step 1: Generate .saf files for site1
        step1="${site1}.saf"
        if ! is_step_complete "$step1"; then
            angsd -bam ${site1}_bamlist.txt -P 64 -minMapQ 30 -minQ 30 -doCounts 1 -GL 1 -doSaf 1 -anc $reference -out ${site1}
            if [ $? -eq 0 ]; then
                log_step "$step1"
            else
                echo "Failed to generate .saf file for $site1"
                exit 1
            fi
        fi

        # Step 2: Generate .saf files for site2
        step2="${site2}.saf"
        if ! is_step_complete "$step2"; then
            angsd -bam ${site2}_bamlist.txt -P 64 -minMapQ 30 -minQ 30 -doCounts 1 -GL 1 -doSaf 1 -anc $reference -out ${site2}
            if [ $? -eq 0 ]; then
                log_step "$step2"
            else
                echo "Failed to generate .saf file for $site2"
                exit 1
            fi
        fi

        # Step 3: Create a 2D SFS prior
        step3="${site1}_${site2}.sfs"
        if ! is_step_complete "$step3"; then
            realSFS ${site1}.saf.idx ${site2}.saf.idx -P 64 > ${site1}_${site2}.sfs
            if [ $? -eq 0 ]; then
                log_step "$step3"
            else
                echo "Failed to create 2D SFS prior for ${site1} and ${site2}"
                exit 1
            fi
        fi

        # Step 4: Generate 2D FST index
        step4="${site1}_${site2}.fst.idx"
        if ! is_step_complete "$step4"; then
            realSFS fst index ${site1}.saf.idx ${site2}.saf.idx -sfs ${site1}_${site2}.sfs -fstout ${site1}_${site2} -P 64
            if [ $? -eq 0 ]; then
                log_step "$step4"
            else
                echo "Failed to generate 2D FST index for ${site1} and ${site2}"
                exit 1
            fi
        fi

        # Step 5: Estimate global FST and save to a file
        step5="${site1}_${site2}_fst.txt"
        if ! is_step_complete "$step5"; then
            realSFS fst stats ${site1}_${site2}.fst.idx -P 64 > ${site1}_${site2}_fst.txt
            if [ $? -eq 0 ]; then
                log_step "$step5"
            else
                echo "Failed to estimate global FST for ${site1} and ${site2}"
                exit 1
            fi
        fi
    done
done
