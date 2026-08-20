#!/bin/bash
#SBATCH -A dewoody
#SBATCH -p highmem
#SBATCH --job-name=hybrid_fast
#SBATCH -N 1
#SBATCH -n 192
#SBATCH --mem=180G
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# build aim genotype inputs for the hybrid triangle analysis (see triangle.R)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
OUT="$PROJ/hybrid_fast"
ENP="$PROJ/lists/east_bams.txt"            # eastern reference (all ENP)
WNPREF="$PROJ/lists/wnp6_extreme_bams.txt" # western reference (low-eastern WNP)
ALL="$PROJ/bams_list.txt"                  # everyone, for plotting
DELTA=0.30                                  # AIM allele-frequency-difference threshold

NCORES=192            # <- raise to 192 if your node/partition allows
T=4                   # threads per ANGSD job
CONC=$(( NCORES / T ))  # concurrent jobs
mkdir -p "$OUT"; cd "$OUT"
export REF T

# 0) candidate sites = genome-wide SNP set (major/minor) from pca.sh; one chr list
zcat "$PROJ/PCA/whale.mafs.gz" | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' > global.sites
angsd sites index global.sites
cut -f1 global.sites | uniq | sort -u > chr_list.txt
echo "chromosomes to scatter over: $(wc -l < chr_list.txt)   cores=$NCORES  T=$T  concurrent=$CONC"

# -------- 1) parental allele frequencies, scattered by chromosome --------
run_maf () {                     # args: pop  chr
  pop="$1"; chr="$2"
  case "$pop" in enp) L="$ENP";; wnp_pure) L="$WNPREF";; esac
  angsd -bam "$L" -ref "$REF" -GL 1 -sites global.sites -r "${chr}:" \
        -doMajorMinor 3 -doMaf 1 -minMapQ 20 -minQ 20 -P "$T" \
        -out "tmp_${pop}_${chr}" > "log_${pop}_${chr}.txt" 2>&1
}
export -f run_maf
export ENP WNPREF

for pop in enp wnp_pure; do
  while read -r chr; do echo "$pop $chr"; done < chr_list.txt
done | xargs -P "$CONC" -n 2 bash -c '
module --force purge
module load biocontainers
module load angsd
run_maf "$@"
' _

# merge per-chromosome mafs (keep a single header)
for pop in enp wnp_pure; do
  first=$(head -1 chr_list.txt)
  zcat "tmp_${pop}_${first}.mafs.gz" | head -1 > "${pop}.mafs"
  while read -r chr; do
    [ -f "tmp_${pop}_${chr}.mafs.gz" ] && zcat "tmp_${pop}_${chr}.mafs.gz" | tail -n +2 >> "${pop}.mafs"
  done < chr_list.txt
  gzip -f "${pop}.mafs"
done

# -------- 2) select AIMs (|Δ freq| >= DELTA) --------
module load anaconda 2>/dev/null || true
python3 - "$DELTA" <<'PY'
import gzip, csv, sys
DELTA=float(sys.argv[1])
def load(f):
    d={}
    with gzip.open(f,'rt') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            d[(row['chromo'],row['position'])]=(row['major'],row['minor'],
                                                float(row['knownEM']),int(float(row['nInd'])))
    return d
enp,wnp=load('enp.mafs.gz'),load('wnp_pure.mafs.gz')
n=0
with open('aims.tsv','w') as out, open('AIM.sites','w') as sites:
    out.write('chr\tpos\tmajor\tminor\tenp_freq\twnp_freq\tdelta\n')
    for k in enp:
        if k not in wnp: continue
        maj,mnr,pe,ne=enp[k]; _,_,pw,nw=wnp[k]
        if ne<16 or nw<4: continue
        d=abs(pe-pw)
        if d>=DELTA:
            out.write(f"{k[0]}\t{k[1]}\t{maj}\t{mnr}\t{pe:.4f}\t{pw:.4f}\t{d:.4f}\n")
            sites.write(f"{k[0]}\t{k[1]}\t{maj}\t{mnr}\n"); n+=1
print(f"AIMs selected (|Δfreq| >= {DELTA}): {n}")
PY
sort -k1,1 -k2,2n AIM.sites -o AIM.sites
angsd sites index AIM.sites
cut -f1 AIM.sites | uniq | sort -u > aim_chr_list.txt

module --force purge
module load biocontainers
module load angsd


# -------- 3) genotype posteriors at AIMs, scattered by chromosome --------
run_geno () {                    # arg: chr
  chr="$1"
  angsd -bam "$ALL" -ref "$REF" -GL 1 -sites AIM.sites -r "${chr}:" \
        -doMajorMinor 3 -doMaf 1 -doPost 2 -doGeno 8 -minMapQ 20 -minQ 20 -P "$T" \
        -out "tmp_geno_${chr}" > "log_geno_${chr}.txt" 2>&1
}
export -f run_geno; export ALL

xargs -P "$CONC" -n 1 bash -c '
module --force purge
module load biocontainers
module load angsd
run_geno "$@"
' _ < aim_chr_list.txt

# concatenate per-chromosome geno posteriors (gzip streams concatenate cleanly;
# triangle.R matches rows to AIMs by chr:pos, so order does not matter)
> tri_geno.geno.gz
while read -r chr; do
  [ -f "tmp_geno_${chr}.geno.gz" ] && cat "tmp_geno_${chr}.geno.gz" >> tri_geno.geno.gz
done < aim_chr_list.txt
sed 's#.*/##; s#\.bam$##' "$ALL" > tri_geno.inds

# tidy scratch
rm -f tmp_*.arg tmp_*.mafs.gz tmp_geno_*.arg tmp_geno_*.mafs.gz

echo "Done. Copy aims.tsv, tri_geno.geno.gz, tri_geno.inds to the R working dir."