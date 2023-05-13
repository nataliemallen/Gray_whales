#!/bin/bash
#SBATCH -A fnrpredator
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

#########################################################################
# this script is identifying repeats, mappability and short scaffolds 
# the following files to use in downstream analyses are created:	  	
# ${REF}.fa (reference file with scaffolds>100kb)							
# ok.bed (regions to analyze in angsd etc)									
# check for repeatmasker file on NCBI to skip that step (*_rm.out.gz )
# https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/
#########################################################################

module load biocontainers
module load genmap/1.3.0
module load bedops
module load bedops/2.4.39
module load bedtools/2.30.0
module load repeatmasker

REF=/scratch/bell/allen715/Gray_whales/Reference/GCA_028021215.1_mEscRob2.pri_genomic.fna

# repeatmasker
# if no repeatmasker file is available run RepeatMasker
RepeatMasker -q -species mammal ${REF}.fa 
# make bed file from NCBI repeatmasker
gunzip GCF_000493695.1_BalAcu1.0_rm.out.gz

cat ${REF}.fa.out|tail -n +4|awk '{print $5,$6,$7,$11}'|sed 's/ /\t/g' \
> repeats.bed

# build an index of the fasta file(s) whose mappability you want to compute
rm -rf index
$GENMAP index -F ${REF}.fa -I index -S 60

# compute mappability
# k = kmer of 100bp
# E = # two mismatches
$GENMAP map -K 100 -E 2 -I index -O mappability -t -w -bg

# sort bed 
$SORTBED -i repeats.bed > repeats_sorted.bed

# make ${REF}.genome
awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ${REF}.fa.fai > ${REF}.genome

# sort genome file
awk '{print $1, $2, $2}' ${REF}.genome > ref2.genome
sed -i 's/ /\t/g' ref2.genome
$SORTBED -i ref2.genome > ref3.genome
awk '{print $1, $2 }' ref3.genome > ${REF}_sorted.genome
sed -i 's/ /\t/g' ${REF}_sorted.genome

# find nonrepeat regions
bedtools complement -i repeats_sorted.bed -g ${REF}_sorted.genome > nonrepeat.bed

# clean mappability file, remove sites with <1 mappability
awk '$4 == 1' mappability.bedgraph > map.bed

awk 'BEGIN {FS="\t"}; {print $1 FS $2 FS $3}' map.bed > mappability.bed

# sort mappability 
$SORTBED -i mappability.bed > mappability2.bed
sed -i 's/ /\t/g' mappability2.bed

# only include sites that are nonrepeats and have mappability ==1
bedtools subtract -a mappability2.bed -b repeats_sorted.bed > map_nonreapeat.bed

# sort file -- by chr then site
bedtools sort -i map_nonreapeat.bed > filter_sorted.bed

# merge overlapping regions
bedtools merge -i filter_sorted.bed > merged.bed

# remove scaffolds shorter than 1Mb
#bioawk -c fastx '{ if(length($seq) > 1000000) { print ">"$name; print $seq }}' \
#${REF}.fa > ${REF}_1Mb.fa
seqtk seq -L 1000000 ${REF}.fa > ${REF}_1Mb.fa

# index
samtools faidx ${REF}_1Mb.fa

# make list with the >1Mb scaffolds
awk '{ print $1, $2, $2 }' ${REF}_1Mb.fa.fai > chrs.info

# replace column 2 with zeros
awk '$2="0"' chrs.info > chrs.bed

# make tab delimited
sed -i 's/ /\t/g' chrs.bed

# make chrs.txt
cut -f 1 chrs.bed > chrs.txt

sort chrs.txt > sorted_chrs.txt

cat satc_sexlinked_scaff.list satc_XZ_scaff.list | uniq | sort -t _ -k 2 -g  > sex.scafs


# remove sex chromosome scaffolds 
comm -1 -3 sex.scafs sorted_chrs.txt > autosomes.txt

xargs samtools faidx ${REF}_1Mb.fa < autosomes.txt > autosomes_1Mb.fa

samtools faidx autosomes_1Mb.fa

# make bed file with the autosomes, 100k, no repeats, mappability =1 sites
awk '{ print $1, $2, $2 }' autosomes_1Mb.fa.fai > autosomes_1Mb.info

# replace column 2 with zeros
awk '$2="0"' autosomes_1Mb.info > autosomes_1Mb.bed

# make tab delimited
sed -i 's/ /\t/g' autosomes_1Mb.bed

# only include scaffolds in merged.bed if they are in autosomes.txt
bedtools intersect -a autosomes_1Mb.bed -b merged.bed > filteredSites.bed	

# get sex chromosome info

xargs samtools faidx ${REF}.fa < idxstats/sex.scafs > satc.fa

samtools faidx satc.fa

# number of sex linked scaffolds
wc -l satc.fa.fai
# length of sex linked scaffolds
cut -f 2 satc.fa.fai | paste -sd+ | bc

# compare sex.scafs and filteredsites
seqtk seq -L 1000000 satc.fa > satc_1Mb.fa
samtools faidx satc_1Mb.fa

cut -f 1 satc_1Mb.fa.fai > sex_1Mb.list

cut -f 1 filteredSitesSorted.bed | uniq | sort -t _ -k 2 -g  > autosome_1Mb.list

grep -wFf autosome_1Mb.list sex_1Mb.list  > removeFromVCF.list

less removeFromVCF.list

# use --not-chr in vcf to filter
# http://vcftools.sourceforge.net/man_latest.html#SITE%20FILTERING%20OPTIONS


# count number of sites in ${REF}_1Mb.fa and ok.bed
# ok.bed == 1373701814 (55.5%)
# ${REF}_1Mb == 2443890318 (98.7%)
# ref == 2476100537

# check that the number of scaffolds in autosomes and ok.bed matches
awk '{ print $1 }' filteredSites.bed > ok.scaffolds
uniq ok.scaffolds > check

awk '{ print $1 }' autosomes_1Mb.bed > auto.scaffolds
uniq auto.scaffolds > check2
wc -l check2
wc -l check

# check if there are scaffolds with <1 mappability for all sites
awk '{ print $1 }' mappability2.bed > mappability.scaffolds
uniq mappability.scaffolds > check3
wc -l check3

# check norepeats file
awk '{ print $1 }' nonrepeat.bed > nonrepeat.scaffolds
uniq nonrepeat.scaffolds > check4
wc -l check4

rm -rf *nonrepeat*
rm -rf chrs.*
rm -rf mappability*bed
rm -rf autosomes*
rm -rf check*
rm -rf ok.*
rm -rf *.genome

# END
