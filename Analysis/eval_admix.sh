#!/bin/bash
#SBATCH -A highmem
#SBATCH --job-name=eval_admix
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 06:000:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

module load biocontainers
module load angsd

# set path to beagle
beagle="/scratch/negishi/allen715/Gray_whales/admixture/whale.beagle.gz"  

/home/allen715/evalAdmix/evalAdmix -beagle "$beagle" -fname admix_K1_run16.fopt.gz -qname admix_K1_run16.qopt -P 32
mv output.corres.txt K1_output.corres.txt

/home/allen715/evalAdmix/evalAdmix -beagle "$beagle" -fname admix_K2_run16.fopt.gz -qname admix_K2_run16.qopt -P 32
mv output.corres.txt K2_output.corres.txt

/home/allen715/evalAdmix/evalAdmix -beagle "$beagle" -fname admix_K3_run16.fopt.gz -qname admix_K3_run16.qopt -P 32
mv output.corres.txt K3_output.corres.txt

/home/allen715/evalAdmix/evalAdmix -beagle "$beagle" -fname admix_K4_run16.fopt.gz -qname admix_K4_run16.qopt -P 32
mv output.corres.txt K4_output.corres.txt

/home/allen715/evalAdmix/evalAdmix -beagle "$beagle" -fname admix_K5_run16.fopt.gz -qname admix_K5_run16.qopt -P 32
mv output.corres.txt K5_output.corres.txt
