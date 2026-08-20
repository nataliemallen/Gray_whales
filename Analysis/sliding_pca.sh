#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=sliding_pca
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 14-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# sliding-window pca in 100 kb windows (Fig S5; plotted in sliding_pca.R)

module --force purge
module load biocontainers
module load angsd
module load pcangsd

PROJ="/scratch/gautschi/allen715/2026_whales"
BEAGLE="$PROJ/PCA/whale.beagle.gz"
OUT="$PROJ/sliding_pca"
WIN=100000          # 100 kb windows
MINSNP=10           # minimum SNPs to attempt a window PCA
mkdir -p "$OUT/windows" "$OUT/cov"; cd "$OUT"

# Save the Beagle header for re-use
zcat "$BEAGLE" | head -n1 > windows/header.txt

# Split the Beagle body into per-window files.
# Marker IDs are CHR_POS where CHR may itself contain underscores, so we take the
# substring after the LAST underscore as the position.
zcat "$BEAGLE" | tail -n +2 | awk -v win="$WIN" -v dir="$OUT/windows" '
{
  m=$1
  p=match(m, /_[0-9]+$/)           # last _<digits>
  chr=substr(m,1,p-1)
  pos=substr(m,p+1)
  w=int(pos/win)*win
  key=chr"@"w
  gsub(/[^A-Za-z0-9._@]/,"_",key)
  print >> (dir"/win_"key".body")
}'

# Run PCAngsd per window with enough SNPs
> window_manifest.txt
for body in windows/win_*.body; do
  n=$(wc -l < "$body")
  [ "$n" -lt "$MINSNP" ] && continue
  base=$(basename "$body" .body)          # win_<chr>@<start>
  tag=${base#win_}
  chr=${tag%@*}; start=${tag#*@}
  cat windows/header.txt "$body" | gzip > "windows/${base}.beagle.gz"
  if pcangsd -b "windows/${base}.beagle.gz" -o "cov/${base}" -t 8 2>/dev/null; then
    echo -e "${chr}\t${start}\t${n}\tcov/${base}.cov" >> window_manifest.txt
  fi
done
echo "Windows with PCA: $(wc -l < window_manifest.txt)"

