###need to modify for whales - run first command as separate script###

#!/bin/sh
#SBATCH -A fnrtowhee
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH --job-name=phylo_tree
#SBATCH -e phylo_tree.err
#SBATCH -o phylo_tree.o

module purge
module load biocontainers
module load vcf-kit/0.2.9
module load bioinfo
module load bcftools/1.15.1
module load iqtree/2.1.2

##WGR tree
#Generate a vcf file for all samples
cd /scratch/bell/eheenken/projects/current_pupfish/final_bams
bcftools mpileup -f /scratch/bell/eheenken/projects/current_pupfish/ref/GCF_016077235.1_ASM1607723v1_genomic.fna \
 -b bam.list | bcftools call -mv -Ov -o /scratch/bell/eheenken/projects/current_pupfish/wgs_tree/vcf_with_var
cd /scratch/bell/eheenken/projects/current_pupfish/wgs_tree
#Generate fasta-alignment from variant calls
vk phylo fasta vcf_with_var
##tree generation -maximum likelihood tree with 1000 boostraps
cd /scratch/bell/eheenken/projects/current_pupfish/wgs_tree
iqtree -s vcf_kit_tree.o -B 1000 -alrt 1000 -T 64
