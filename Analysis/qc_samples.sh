#!/bin/bash
# qc step 1: per-sample checks (metadata, alignment, heterozygosity, within-bam mixture)
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 7:00:00
#SBATCH --job-name=qc1_persample
#SBATCH -e logs/%x_%A_%a.err
#SBATCH -o logs/%x_%A_%a.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

set -uo pipefail
module --force purge
module load biocontainers
module load samtools
module load bcftools
mkdir -p logs

# config
bam_list="/scratch/gautschi/allen715/whale/bams_list.txt"      # one BAM path per line
ref_fasta="/scratch/gautschi/allen715/2026_whales/reference/ref.fa"
outdir="/scratch/gautschi/allen715/2026_whales/QC"
THREADS=${SLURM_CPUS_PER_TASK:-16}
MINMAPQ=30
MINQ=20

mkdir -p "$outdir/persample" "$outdir/het" "$outdir/mixture" "$outdir/tmp"
[ -f "${ref_fasta}.fai" ] || samtools faidx "$ref_fasta"
GENOME=$(awk '{s+=$2} END{print s}' "${ref_fasta}.fai")

# pick this task's BAM
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$bam_list")
[ -z "$bam" ] && { echo "no bam at line ${SLURM_ARRAY_TASK_ID}"; exit 1; }
name=$(basename "$bam" .bam)
echo ">>> [$SLURM_ARRAY_TASK_ID] $name"
[ -f "${bam}.bai" ] || [ -f "${bam%.bam}.bai" ] || samtools index -@ $THREADS "$bam"

# (1) metadata audit
hdr="$outdir/persample/$name.header.txt"
samtools view -H "$bam" > "$hdr"
n_rg=$(grep -c '^@RG' "$hdr" || true)
n_bwa=$(grep -c '^@PG.*ID:bwa' "$hdr" || true)
distinct_sm=$(grep '^@RG' "$hdr" | grep -o 'SM:[^[:space:]]*' | sort -u | sed 's/SM://' | paste -sd';' -)
n_sm=$(grep '^@RG' "$hdr" | grep -o 'SM:[^[:space:]]*' | sort -u | wc -l)
n_lb=$(grep '^@RG' "$hdr" | grep -o 'LB:[^[:space:]]*' | sort -u | wc -l)

# (2) alignment qc
stats="$outdir/persample/$name.stats.txt"
samtools stats -@ $THREADS "$bam" > "$stats"
error_rate=$(grep '^SN' "$stats" | grep 'error rate:'                    | cut -f3)
ins_size=$(  grep '^SN' "$stats" | grep 'insert size average:'           | cut -f3)
ins_sd=$(    grep '^SN' "$stats" | grep 'insert size standard deviation:'| cut -f3)
read_len=$(  grep '^SN' "$stats" | grep 'average length:'                | cut -f3)
raw=$(       grep '^SN' "$stats" | grep 'raw total sequences:'           | cut -f3)
dup=$(       grep '^SN' "$stats" | grep 'reads duplicated:'              | cut -f3)
fs=$(samtools flagstat -@ $THREADS "$bam")
pct_mapped=$(echo "$fs" | awk '/ mapped \(/{gsub(/[()%]/,"",$6);print $6;exit}')
pct_pp=$(    echo "$fs" | awk '/properly paired \(/{gsub(/[()%]/,"",$6);print $6;exit}')
pct_dup=$(awk -v d="$dup" -v r="$raw" 'BEGIN{if(r>0)printf "%.3f",100*d/r; else print "NA"}')

# depth + breadth in one pass (all positions)
read mean_depth breadth1 breadth4 < <(
  samtools depth -a -q $MINQ -Q $MINMAPQ "$bam" \
  | awk -v G="$GENOME" '{s+=$3; if($3>=1)c1++; if($3>=4)c4++}
        END{ if(G>0) printf "%.3f %.4f %.4f", s/G, c1/G, c4/G; else print "NA NA NA" }'
)

# (3) heterozygosity (genome-wide folded sfs)
angsd -i "$bam" -ref "$ref_fasta" -anc "$ref_fasta" \
      -GL 1 -doSaf 1 -minMapQ $MINMAPQ -minQ $MINQ -C 50 -baq 1 \
      -nThreads $THREADS -out "$outdir/het/$name"
realSFS "$outdir/het/$name.saf.idx" -fold 1 -P $THREADS > "$outdir/het/$name.sfs" 2>/dev/null
het=$(awk '{t=$1+$2+$3; if(t>0)printf "%.6f",$2/t; else print "NA"}' "$outdir/het/$name.sfs")

# (4) within-bam mixture (split by read group -> ibs)
mix_med="NA"; mix_max="NA"; n_units=$n_rg
if [ "${n_rg:-0}" -gt 1 ]; then
  ud="$outdir/mixture/$name"; mkdir -p "$ud/units"
  samtools split -@ $THREADS -f "$ud/units/%!.bam" "$bam" 2>/dev/null || true
  ls "$ud"/units/*.bam > "$ud/units.list" 2>/dev/null || true
  nU=$(wc -l < "$ud/units.list" 2>/dev/null || echo 0)
  n_units=$nU
  if [ "$nU" -gt 1 ]; then
    for u in "$ud"/units/*.bam; do samtools index "$u"; done
    angsd -bam "$ud/units.list" -ref "$ref_fasta" \
          -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
          -doIBS 1 -makeMatrix 1 -doCounts 1 \
          -minMapQ $MINMAPQ -minQ $MINQ -nThreads $THREADS -out "$ud/ibs" 2>/dev/null
    if [ -f "$ud/ibs.ibsMat" ]; then
      read mix_med mix_max < <(
        awk '{for(j=1;j<=NF;j++) if(j!=NR){v[++k]=$j}}
             END{ if(k==0){print "NA NA";exit}
                  asort(v); med=(k%2)?v[(k+1)/2]:(v[k/2]+v[k/2+1])/2;
                  printf "%.4f %.4f", med, v[k] }' "$ud/ibs.ibsMat"
      )
    fi
  fi
fi

# write one row
row="$outdir/persample/$name.qc.tsv"
echo -e "sample\tn_RG\tn_bwa_runs\tn_distinct_SM\tdistinct_SM\tn_distinct_LB\tn_units\tmean_depth\tbreadth_1x\tbreadth_4x\tpct_dup\terror_rate\tins_size\tins_sd\tread_len\tpct_mapped\tpct_proper_pair\theterozygosity\tmix_ibs_median\tmix_ibs_max" > "$row"
echo -e "${name}\t${n_rg}\t${n_bwa}\t${n_sm}\t${distinct_sm}\t${n_lb}\t${n_units}\t${mean_depth}\t${breadth1}\t${breadth4}\t${pct_dup}\t${error_rate}\t${ins_size}\t${ins_sd}\t${read_len}\t${pct_mapped}\t${pct_pp}\t${het}\t${mix_med}\t${mix_max}" >> "$row"

echo ">>> done $name"

