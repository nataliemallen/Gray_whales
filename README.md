# Gray_whales

# Mapping 
Begin with folder of FASTQ files
Create sample.list file in reads directory (lists 1 file name for each pair of reads, format files so that only ends differ with R1/R2 - check example_sample.list)
Use parse_genome.sh, reference_format.sh, and dictionary.sh to prepare reference genome
Use align_whales.sh to create slurm jobs for mapping
Run line in run_jobs.txt to map
Use merge_whales.slurm to merge files from multiple lanes
