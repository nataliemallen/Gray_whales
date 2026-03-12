#!/bin/bash
#SBATCH --job-name=cnv_long
#SBATCH -A fnrdewoody
#SBATCH -t 12-00:00:00
#SBATCH -p cpu
#SBATCH -n 128

REFERENCE_GENOME="/scratch/negishi/allen715/Gray_whales/new_ref/ref.fa"
BAM_DIR="/scratch/negishi/allen715/Gray_whales/final_bams/merged/final_merged_use/"
OUTPUT_DIR="/scratch/negishi/allen715/Gray_whales/CNVs/"
POPULATION1_LIST="east_bams.txt"
POPULATION2_LIST="west_bams_33.txt"

module load biocontainers
module load samtools
module load anaconda
export PATH="/home/allen715/CNVpytor/cnvpytor:$PATH"

mkdir -p ${OUTPUT_DIR}/pytor_files
mkdir -p ${OUTPUT_DIR}/results

process_bam () {
    ID=$1
    BAM_FILE="${BAM_DIR}/${ID}.bam"
    PYTOR_FILE="${OUTPUT_DIR}/pytor_files/${ID}.pytor"

    [ -f "${BAM_FILE}.bai" ] || samtools index "$BAM_FILE"

    cnvpytor -root "${PYTOR_FILE}" -rd "$BAM_FILE"
    cnvpytor -root "${PYTOR_FILE}" -gc "$REFERENCE_GENOME"
    cnvpytor -root "${PYTOR_FILE}" -his 10000 100000
    cnvpytor -root "${PYTOR_FILE}" -partition 10000 100000

    echo "call 10000"  | cnvpytor -root "${PYTOR_FILE}" > "${OUTPUT_DIR}/results/${ID}_calls.10000.tsv"
    echo "call 100000" | cnvpytor -root "${PYTOR_FILE}" > "${OUTPUT_DIR}/results/${ID}_calls.100000.tsv"
}

export -f process_bam
export BAM_DIR OUTPUT_DIR REFERENCE_GENOME

cat $POPULATION1_LIST | xargs -I{} bash -c 'process_bam "$@"' _ {}
cat $POPULATION2_LIST | xargs -I{} bash -c 'process_bam "$@"' _ {}
