#!/bin/bash
# identify and drop the contaminant read-group unit in 038046 by single-read ibs,
# rebuild a cleaned bam from the kept read groups, then revalidate with ngsrelate
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 5-00:00:00
#SBATCH --job-name=clean038046
#SBATCH -e logs/%x_%j.err
#SBATCH -o logs/%x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

set -uo pipefail
module --force purge
module load biocontainers
module load angsd
module load pcangsd
module load samtools
module load bcftools

module load anaconda
conda activate ngsrelate_env 

# module load samtools angsd
mkdir -p logs

# config
target_bam="/scratch/gautschi/allen715/whale/final_bams/merged/final_merged_use/038046_ER-17-0237.bam"
bam_list="/scratch/gautschi/allen715/whale/bams_list.txt"   # all 74, one per line
ref_fasta="/scratch/gautschi/allen715/2026_whales/reference/ref.fa"
outdir="/scratch/gautschi/allen715/2026_whales/QC"
pop_csv="/scratch/gautschi/allen715/2026_whales/whale_nobo_pop.csv"   # ind,sample_ID,pop (same order as bam_list!)

THREADS=${SLURM_CPUS_PER_TASK:-96}
MINMAPQ=30
MINQ=20
DUP_FRAC=0.40     # IBS < DUP_FRAC * (unrelated median) => duplicate-level match

mkdir -p "$outdir/units"
[ -f "${ref_fasta}.fai" ] || samtools faidx "$ref_fasta"
tname=$(basename "$target_bam" .bam)

# 1) split target into per-read-group units
echo ">>> splitting $tname by read group"
samtools split -@ $THREADS -f "$outdir/units/%!.bam" "$target_bam"
for u in "$outdir"/units/*.bam; do samtools index "$u"; done

echo -e "unit_RG_ID\tmapped_reads\tmean_depth" > "$outdir/unit_stats.tsv"
for u in "$outdir"/units/*.bam; do
  rg=$(basename "$u" .bam)
  mr=$(samtools view -c -F 0x904 "$u")
  md=$(samtools depth -a -q $MINQ -Q $MINMAPQ "$u" | awk '{s+=$3;n++} END{if(n)printf "%.3f",s/n; else print 0}')
  echo -e "${rg}\t${mr}\t${md}" >> "$outdir/unit_stats.tsv"
done
echo "---- per-unit stats ----"; column -t "$outdir/unit_stats.tsv"

# 2) ibs matrix: 7 units + all other dataset samples
ibslist="$outdir/ibs_bams.txt"
ls "$outdir"/units/*.bam > "$ibslist"                          # units first
grep -v "$(basename "$target_bam")" "$bam_list" >> "$ibslist"  # then the 73 others
nUNITS=$(ls "$outdir"/units/*.bam | wc -l)
echo ">>> IBS over $(wc -l < "$ibslist") samples ($nUNITS units + others)"

angsd -bam "$ibslist" -ref "$ref_fasta" -GL 1 \
      -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
      -doIBS 1 -makeMatrix 1 -doCounts 1 \
      -minMapQ $MINMAPQ -minQ $MINQ -nThreads $THREADS -out "$outdir/ibs"

# 3) decide which units to keep / drop
python3 - "$outdir/ibs.ibsMat" "$ibslist" "$nUNITS" "$DUP_FRAC" "$outdir" <<'PY'
import sys
mat,lst,nU,frac,outdir=sys.argv[1],sys.argv[2],int(sys.argv[3]),float(sys.argv[4]),sys.argv[5]
names=[l.strip() for l in open(lst) if l.strip()]
base=[n.split('/')[-1][:-4] for n in names]
M=[[float(x) for x in line.split()] for line in open(mat) if line.strip()]
units=list(range(nU)); others=list(range(nU,len(base)))
# reference "unrelated" level = median ibs among the other dataset samples
off=[M[i][j] for i in others for j in others if i<j]
off.sort(); ref=off[len(off)//2]
dupthr=frac*ref
print(f"\nunrelated-median IBS = {ref:.4f}   duplicate threshold = {dupthr:.4f}\n")
keep,drop=[],[]
print("unit decisions:")
for u in units:
    # nearest other dataset sample to this unit
    nn=min(others,key=lambda j:M[u][j]); d=M[u][nn]
    # mean distance to the other units
    du=[M[u][v] for v in units if v!=u]; mu=sum(du)/len(du) if du else float('nan')
    if d < dupthr:
        drop.append(u)
        print(f"  DROP {base[u]:18s} nearest sample={base[nn]} IBS={d:.4f}  (duplicate-level -> CONTAMINANT, source likely {base[nn]})")
    else:
        keep.append(u)
        print(f"  KEEP {base[u]:18s} nearest sample={base[nn]} IBS={d:.4f}  (background-level -> true 038046)")
print(f"\nmean IBS among units: {sum(M[i][j] for i in units for j in units if i<j)/max(1,len(units)*(len(units)-1)//2):.4f}")
open(f"{outdir}/keep_rg.txt","w").write("\n".join(base[u] for u in keep)+"\n")
open(f"{outdir}/drop_rg.txt","w").write("\n".join(base[u] for u in drop)+"\n")
if not drop:
    print("\n*** NO unit is a duplicate of another sample. Contamination is NOT")
    print("    separable by read group -> DROP the whole sample, do not 'clean'. ***")
elif not keep:
    print("\n*** ALL units look like duplicates of other samples -> drop whole sample. ***")
else:
    print(f"\n-> keeping {len(keep)} read groups, dropping {len(drop)}. See keep_rg.txt / drop_rg.txt")
PY

# 4) rebuild cleaned bam from kept read groups (no re-alignment)
if [ -s "$outdir/drop_rg.txt" ] && [ -s "$outdir/keep_rg.txt" ]; then
  clean="$outdir/${tname}.cleaned.bam"
  samtools view -h -@ $THREADS -R "$outdir/keep_rg.txt" -b -o "$clean" "$target_bam"
  samtools index "$clean"
  echo ">>> cleaned BAM: $clean"
  echo "    reads before: $(samtools view -c -F 0x904 "$target_bam")  after: $(samtools view -c -F 0x904 "$clean")"

  # 5) validate: re-run ngsrelate with cleaned bam swapped in
  vlist="$outdir/validate_bams.txt"
  sed "s|$target_bam|$clean|" "$bam_list" > "$vlist"
  N=$(wc -l < "$vlist")
  cd "$outdir"
  angsd -bam "$vlist" -ref "$ref_fasta" -GL 2 -doMajorMinor 1 -doMaf 1 \
        -SNP_pval 1e-6 -minMaf 0.05 -doGlf 3 -minMapQ $MINMAPQ -minQ $MINQ \
        -nThreads $THREADS -out validate_glf3
  zcat validate_glf3.mafs.gz | cut -f6 | sed 1d > validate_freq
  "$ngsrelate" -g validate_glf3.glf.gz -n "$N" -f validate_freq -O validate.res -p $THREADS

  python3 - validate.res "$vlist" "$clean" <<'PY'
import sys,csv
res,lst,clean=sys.argv[1],sys.argv[2],sys.argv[3]
names=[l.strip().split('/')[-1][:-4] for l in open(lst) if l.strip()]
ti=names.index(clean.split('/')[-1][:-4])
import statistics as st
R0=[];K=[];nrel=0
with open(res) as fh:
    first=fh.readline(); d='\t' if '\t' in first else ','; fh.seek(0)
    rd=csv.reader(fh,delimiter=d); h=[x.strip() for x in next(rd)]
    a=h.index('a') if 'a' in h else 0; b=h.index('b') if 'b' in h else 1
    iR0=h.index('R0'); iK=h.index('KING'); iT=h.index('theta')
    for p in rd:
        try: ai=int(p[a]); bi=int(p[b])
        except: continue
        if ti in (ai,bi):
            R0.append(float(p[iR0])); k=float(p[iK]); t=float(p[iT]); K.append(k)
            if k>=0.0884 and t<0.0442: nrel+=1
print("\n===== VALIDATION (cleaned 038046) =====")
print(f"median R0 = {st.median(R0):.3f}  (contaminated was ~0.19; healthy ~0.55-0.60)")
print(f"median KING = {st.median(K):.3f}  (contaminated was ~0.12; healthy ~ -0.02)")
print(f"'related-to-many' partners = {nrel}  (contaminated was 67; healthy ~0-1)")
print("PASS if R0 back near 0.55-0.60 and related-partners ~0." )
PY
else
  echo ">>> No cleaning performed (see decision above)."
fi

