### from Andrew Black
### keeps only chromosomes in bam files, removes mitochondrial genomes and unmapped scaffolds

#!/bin/bash
#SBATCH -A fnrchook
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --error=align.err
#SBATCH --output=align.out
#SBATCH --job-name=produce_align_SLURMM_jobs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=



#Move to fastq containing directory
cd /scratch/bell/dewoody/gray_whale_final_bams_allen
mkdir jobs2
#Define variables to shorten commands

while read -a line
do 
        echo "#!/bin/sh -l
#SBATCH -A fnrchook
#SBATCH -t 05-00:00:00
#SBATCH --job-name=${line[0]}_aln
#SBATCH --error=${line[0]}_aln.e
#SBATCH --output=${line[0]}_aln.o
#SBATCH --mem=30G
module --force purge
module load bioinfo
module load samtools
#Move to the bam containing folder
cd /scratch/bell/dewoody/gray_whale_final_bams_allen


#samtools index  ${line[0]}
#samtools depth -a ${line[0]} \
#| awk '{c++;s+=\$3}END{print s/c}' \
#> ${line}.post.meandepth.txt

#samtools depth -a ${line[0]}  \
#| awk '{c++; if(\$3>0) total+=1}END{print (total/c)*100}' \
#> ${line}.post.1xbreadth.txt

samtools view -L chrom.bed -h ${line[0]} > ${line[0]}_chrom.bam
echo done" > ./jobs2/${line[0]}_aln.sh

done < ./depth


#Print out contents of depth and breadth
#for i in `cat depth` ; do cat $i.post.1xbreadth.txt ; done 
#for i in `cat depth` ; do cat $i.post.meandepth.txt ; done 
