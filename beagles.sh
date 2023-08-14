### from Andrew Black

#!/bin/bash
#SBATCH -A fnrchook
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --error=beagle.err
#SBATCH --output=beagle.out
#SBATCH --job-name=produce_regions_SLURMM_jobs


mkdir jobs_beagle
while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A fnrchook
#SBATCH -n 10
#SBATCH -t 10-00:00:00
#SBATCH --job-name=Beagle_${line[0]}
#SBATCH --error=Beagle_${line[0]}.e
#SBATCH --output=Beagle_${line[0]}.o
#SBATCH --mem=20G
module load biocontainers
module load angsd

#Move to bam folder
mkdir /scratch/bell/dewoody/gray_whale_final_bams_allen/BEAGLE/
cd /scratch/bell/dewoody/gray_whale_final_bams_allen
angsd -GL 1 -out BEAGLE/${line[0]} -r ${line[0]} -minQ 30 -P 10 \
-minInd 53 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -skipTriallelic 1 -SNP_pval 1e-6 -bam bam.list  \
-ref /scratch/bell/dewoody/gray_whale_final_bams_allen/original.fa "  > ./jobs_beagle/filtered_${line[0]}.sh
done < ./regions.txt

#for i in `ls -1 *sh`; do  echo "sbatch $i" ; done > jobs ; source ./jobs

#cd /scratch/bell/dewoody/gray_whale_final_bams_allen/BEAGLE/
#for f in *.beagle.gz; do zcat ${f} | sed 1d   > tmpfile; mv tmpfile ${f}_trim; done
#zcat NW_026294820.1.beagle.gz | head -n 1 > beagle.header
#cat beagle.header  *_trim > filtered.beagle ; gzip filtered.beagle 
