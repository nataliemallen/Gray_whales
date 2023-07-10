#!/bin/bash
#SBATCH --job-name=align_whales
#SBATCH -A fnrpredator
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12-00:00:00
#SBATCH --error=align.err
#SBATCH --output=align.out
#SBATCH --job-name=produce_align_SLURMM_jobs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools/2.18.2
module load bedops
module load GATK/3.6.0
module load samtools/1.5

#Make sample list

cd /scratch/bell/allen715/Gray_whales/Reads/
#Make directory to hold all SLURMM jobs
mkdir jobs

#Define variables to shorten commands
REF=/scratch/bell/allen715/Gray_whales/Reference/original.fa
DICT=/scratch/bell/allen715/Gray_whales/Reference/original.fa.dict

#bwa index $REF
#samtools faidx $REF
#PicardCommandLine CreateSequenceDictionary reference=$REF output=$DICT


while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A fnrpredator
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 05-00:00:00
#SBATCH --job-name=${line[0]}_align_stats
#SBATCH --error=${line[0]}_align_stats.e
#SBATCH --output=${line[0]}_align_stats.o
#SBATCH --mem=20G

module --force purge
module load bioinfo
module load bwa
module load picard-tools/2.18.2
module load bedops
module load GATK/3.6.0
module load samtools/1.5


#Move to the paired-end fastq containing folder
cd  /scratch/bell/allen715/Gray_whales/Reads/

# Align sample to indexed reference genome
bwa mem -t 10 -M -R \"@RG\tID:group1\tSM:${line[0]}\tPL:illumina\tLB:lib1\tPU:unit1\" \
/scratch/bell/allen715/Gray_whales/Reference/original.fa \
${line[0]}_R1_001.fastq.gz ${line[0]}_R2_001.fastq.gz > ../aligned/${line[0]}.sam

#Move to aligned directory
cd ../aligned

#Validate sam file
PicardCommandLine ValidateSamFile I=${line[0]}.sam MODE=SUMMARY O=${line[0]}.sam.txt

#Sort validated sam file by read coordinate
PicardCommandLine SortSam INPUT=${line[0]}.sam OUTPUT=sorted_${line[0]}.bam SORT_ORDER=coordinate

#Get summary stats on initial alignments:
samtools flagstat sorted_${line[0]}.bam > ${line[0]}_mapped.txt

#Mark duplicates
PicardCommandLine MarkDuplicates INPUT=sorted_${line[0]}.bam OUTPUT=dedup_${line[0]}.bam METRICS_FILE=metrics_${line[0]}.bam.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500

#Index in prep for realignment
PicardCommandLine BuildBamIndex INPUT=dedup_${line[0]}.bam

# local realignment of reads
GenomeAnalysisTK -T RealignerTargetCreator -nt 10 -R /scratch/bell/allen715/Gray_whales/Reference/original.fa -I dedup_${line[0]}.bam -o ${line[0]}_forIndelRealigner.intervals

#Realign with established intervals
GenomeAnalysisTK -T IndelRealigner -R /scratch/bell/allen715/Gray_whales/Reference/original.fa -I dedup_${line[0]}.bam -targetIntervals ${line[0]}_forIndelRealigner.intervals -o ${line[0]}_indel.bam

#Make new directory
mkdir ../final_bams

#Fix mate info
PicardCommandLine FixMateInformation INPUT=dedup_${line[0]}.bam OUTPUT=${line[0]}.fixmate.bam SO=coordinate CREATE_INDEX=true

#   Remove unmapped (4), secondary (256), QC failed (512), duplicate (1024), and
#   supplementary (2048) reads from indel-realigned BAMs, and keep only reads
#   mapped in a proper pair (2) to regions in a BED file (produced from QC_reference.sh)

samtools view -@ 10 -q 30 -b -F 3844 -f 2 -L /scratch/bell/allen715/Gray_whales/Reference/ok.bed ${line[0]}.fixmate.bam > ../final_bams/${line[0]}_filt.bam 

#Move into the final directory
cd ../final_bams/
#Index bam file
PicardCommandLine BuildBamIndex INPUT=${line[0]}_filt.bam

#Summary stats on files
samtools depth -a ${line[0]}_filt.bam | awk '{c++;s+=$3}END{print s/c}' > ${line[0]}_depth.txt
samtools depth -a ${line[0]}_filt.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${line[0]}_1x_breadth.txt

samtools stats ${line[0]}_filt.bam > ${line[0]}_samtools_stats.txt" > ./jobs/${line[0]}_alignment.sh

done < ./sample.list

#when finished, paste below into jobs directory to run all slurm jobs:
#for i in `ls -1 *sh`; do  echo "sbatch $i" ; done > slurmm_jobs ; source ./slurmm_jobs
