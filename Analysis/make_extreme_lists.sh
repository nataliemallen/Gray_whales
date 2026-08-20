#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=fst_extreme
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 4-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# build the extreme ENP/WNP bam lists (6 each) for fst_extreme.sh

PROJ="/scratch/gautschi/allen715/2026_whales"
BAMDIR="$PROJ/bams"
LISTS="$PROJ/lists"
SEED=42
mkdir -p "$LISTS"

# --- WNP: the 6 chosen "extreme" west individuals ---
WNP=(
  038001_ER-14-0163
  038007_ER-16-0047
  038020_ER-16-0063
  038023_ER-16-0067
  038026_ER-16-0073
  038057_ER-14-0165
)
: > "$LISTS/wnp6_extreme_bams.txt"
for s in "${WNP[@]}"; do
  printf '%s\n' "$BAMDIR/$s.bam" >> "$LISTS/wnp6_extreme_bams.txt"
done

# --- ENP: 6 random east individuals (from lists/east.txt, seeded) ---
python3 - "$LISTS/east.txt" "$SEED" "$BAMDIR" "$LISTS/enp6_bams.txt" <<'PY'
import random, sys
src, seed, bamdir, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
ids = [l.strip() for l in open(src) if l.strip()]
random.seed(seed)
pick = sorted(random.sample(ids, 6))
with open(out, "w") as fh:
    for s in pick:
        fh.write(f"{bamdir}/{s}.bam\n")
PY

echo "WNP (6, fixed):";   cat "$LISTS/wnp6_extreme_bams.txt"
echo "ENP (6, seed=$SEED):"; cat "$LISTS/enp6_bams.txt"