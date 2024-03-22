#!/bin/bash
#SBATCH --job-name=heterozygosity
#SBATCH -A highmem
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

mkdir jobs_het
mkdir HET

while read -a line
do 
	echo "#!/bin/bash
#SBATCH -A highmem
#SBATCH -n 10
#SBATCH -t 1-00:00:00
#SBATCH --job-name=${line[0]}_het_stats
#SBATCH --error=${line[0]}_het_stats.e
#SBATCH --output=${line[0]}_het_stats.o

module load biocontainers
module load angsd

#Move to the bams folder
cd /scratch/negishi/allen715/Gray_whales/final_bams/merged/

angsd -i ${line[0]}.bam -ref /scratch/negishi/allen715/Gray_whales/reference/ref.fa  -anc /scratch/negishi/allen715/Gray_whales/reference/ref.fa  -dosaf 1 -minMapQ 30 -GL 1 -P 26 -out /scratch/negishi/allen715/Gray_whales/heterozygosity/HET/${line[0]} -doCounts 1 -setMinDepth 3

realSFS -P 10 -fold 1 /scratch/negishi/allen715/Gray_whales/heterozygosity/HET/${line[0]}.saf.idx > /scratch/negishi/allen715/Gray_whales/heterozygosity/HET/${line[0]}_est.ml" > ./jobs_het/${line[0]}_alignment.sh

done < ./sample.list

#for i in `ls -1 *sh`; do  echo "sbatch $i" ; done > jobs ; source ./jobs

#Get individual heterozygosity, proportion of heterozygotes:
#cat ./*ml

#Use output from cat command for calculate prop heterozygote
#DONE
