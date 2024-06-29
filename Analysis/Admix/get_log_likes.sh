#!/bin/bash
#SBATCH -A fnrdewoody
#SBATCH --job-name=admix
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# Directory containing the log files
log_dir="//scratch/negishi/allen715/Gray_whales/admixture/ADX/"

# Output file
output_file="log_likelihoods.txt"

# Empty the output file if it exists
#> $output_file

# Loop through all log files in the directory
for log_file in "$log_dir"/*.log; do
  # Extract the run name from the file name
  run_name=$(grep "outfiles=" "$log_file" | sed 's/.*outfiles=//')
  
  # Extract the log likelihood value
  log_likelihood=$(grep "best like=" "$log_file" | sed 's/.*best like=//;s/ after.*//')
  
  # Write the run name and log likelihood value to the output file
  echo "$run_name $log_likelihood" >> $output_file
done
