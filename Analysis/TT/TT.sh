#!/bin/bash
#SBATCH --job-name=TTbatch
#SBATCH -A fnrdewoody
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=allen715@purdue.edu

module load biocontainers
module load bcftools
module load htslib
module load r

WORKDIR="/scratch/negishi/allen715/Gray_whales/TT-method"
VCF="$WORKDIR/TTo.vcf"
OUTMAIN="$WORKDIR"
SITES="$OUTMAIN/sites/ancestral_supported_sites.txt.bgz"
RSCRIPT_PATH="$WORKDIR/tt_calc.R"
PAIR_FILE="$WORKDIR/TT_pairs.csv"  
RESULTS="$WORKDIR/TT_reslts.csv"

MU=1.11e-8
G=18

# header for the output
echo "Individual1,Individual2,T1,T2,Na,mu (assumed),g (assumed),alpha1,alpha2,theta,t1,t2,v1,v2" > "$RESULTS"

# read CSV (excluding header)
tail -n +2 "$PAIR_FILE" | while IFS=, read -r S1 S2; do
  PAIR="${S1}_${S2}"
  echo "Processing pair: $PAIR"

  # 2dsfs
  OUTPUT_DIR="$WORKDIR/2dsfs"
  SFS_OUTPUT="$OUTPUT_DIR/${PAIR}_unfolded.2dsfs"
  mkdir -p "$OUTPUT_DIR"

  bcftools view -s "$S1,$S2" -T "$SITES" "$VCF" -o temp_samples.vcf
  bcftools query -f '[%GT\n]' -s "$S1" temp_samples.vcf > temp_${S1}_genotypes.txt
  bcftools query -f '[%GT\n]' -s "$S2" temp_samples.vcf > temp_${S2}_genotypes.txt
  paste temp_${S1}_genotypes.txt temp_${S2}_genotypes.txt > temp_genotypes.txt

  awk '
  {
    if ($1 == "./." || $2 == "./.") next;
    geno_combo = $1"_"$2;
    count[geno_combo]++;
  }
  END {
    for (combo in count) {
      print combo, count[combo];
    }
  }' temp_genotypes.txt > "$SFS_OUTPUT"

  rm temp_samples.vcf temp_${S1}_genotypes.txt temp_${S2}_genotypes.txt temp_genotypes.txt

  # run R script
  OUTPREFIX="$WORKDIR/estimates/tt/${PAIR}"
  mkdir -p "$WORKDIR/estimates/tt"
  Rscript "$RSCRIPT_PATH" "$SFS_OUTPUT" "$OUTPREFIX" "$MU" "$G"

  # parse output files
  SPLIT_FILE="${OUTPREFIX}.tt.split.years"
  PARAM_FILE="${OUTPREFIX}.tt.params.res"

  if [[ -f "$SPLIT_FILE" && -f "$PARAM_FILE" ]]; then
    T1=$(grep '^T1' "$SPLIT_FILE" | awk '{print $2}')
    T2=$(grep '^T2' "$SPLIT_FILE" | awk '{print $2}')
    Na=$(grep '^Na' "$SPLIT_FILE" | awk '{print $2}')
    mu=$(grep '^mu' "$SPLIT_FILE" | awk '{print $2}')
    g=$(grep '^g' "$SPLIT_FILE" | awk '{print $2}')

    alpha1=$(grep '^alpha1' "$PARAM_FILE" | awk '{print $2}')
    alpha2=$(grep '^alpha2' "$PARAM_FILE" | awk '{print $2}')
    theta=$(grep '^theta' "$PARAM_FILE" | awk '{print $2}')
    t1=$(grep '^t1' "$PARAM_FILE" | awk '{print $2}')
    t2=$(grep '^t2' "$PARAM_FILE" | awk '{print $2}')
    v1=$(grep '^v1' "$PARAM_FILE" | awk '{print $2}')
    v2=$(grep '^v2' "$PARAM_FILE" | awk '{print $2}')

    echo "$S1,$S2,$T1,$T2,$Na,$mu,$g,$alpha1,$alpha2,$theta,$t1,$t2,$v1,$v2" >> "$RESULTS"
  else
    echo "$S1,$S2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "$RESULTS"
    echo "Warning: Missing result files for $PAIR"
  fi

done

