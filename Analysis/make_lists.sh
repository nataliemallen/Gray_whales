#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=aggregate_qc
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL

# build cohort, east/west, and reproducible downsampled sample lists

set -euo pipefail

PROJ="/scratch/gautschi/allen715/2026_whales"
BAMDIR="$PROJ/bams"
LISTS="$PROJ/lists"
META="$LISTS/sample_pop.tsv"          # cols: sample  pop  breadth  depth
SEED=42                               # fixed seed -> reproducible downsample

# Samples to exclude
declare -A EXCLUDE=(
    ["038046_ER-17-0237"]=1
    ["038003_ER-14-0171"]=1
    ["038016_ER-16-0058"]=1
)

mkdir -p "$LISTS"

# --- population lookup from metadata (exact string keys, no numeric coercion) ---
declare -A POP
while IFS=$'\t' read -r sample pop _; do
  [ "$sample" = "sample" ] && continue   # skip header
  [ -z "$sample" ] && continue
  POP["$sample"]="$pop"
done < "$META"

# --- start clean ---
: > "$PROJ/bams_list.txt"
: > "$LISTS/east.txt"
: > "$LISTS/west.txt"
: > "$LISTS/east_bams.txt"
: > "$LISTS/west_bams.txt"

# --- walk the BAM directory (ground truth for names) ---
shopt -s nullglob
for bam in "$BAMDIR"/*.bam; do
  sample="$(basename "$bam" .bam)"

  # Skip excluded samples
  [[ -n "${EXCLUDE[$sample]:-}" ]] && continue

  printf '%s\n' "$bam" >> "$PROJ/bams_list.txt"

  case "${POP[$sample]:-}" in
    east)
      printf '%s\n' "$sample" >> "$LISTS/east.txt"
      printf '%s\n' "$bam" >> "$LISTS/east_bams.txt"
      ;;
    west)
      printf '%s\n' "$sample" >> "$LISTS/west.txt"
      printf '%s\n' "$bam" >> "$LISTS/west_bams.txt"
      ;;
    *)
      echo "WARNING: no population for '$sample' (not in $META)" >&2
      ;;
  esac
done
shopt -u nullglob

NE=$(wc -l < "$LISTS/east.txt")
NW=$(wc -l < "$LISTS/west.txt")

# --- reproducible downsample of west to match east (n = NE) ---
python3 - "$LISTS/west.txt" "$NE" "$SEED" "$LISTS/west_ds.txt" <<'PY'
import random, sys

src, n, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
ids = [l.rstrip("\n") for l in open(src) if l.strip()]
random.seed(seed)
keep = sorted(random.sample(ids, n))
open(out, "w").write("\n".join(keep) + "\n")
PY

# --- downsampled west BAM paths ---
: > "$LISTS/west_ds_bams.txt"
while IFS= read -r s; do
  [ -z "$s" ] && continue
  printf '%s\n' "$BAMDIR/$s.bam" >> "$LISTS/west_ds_bams.txt"
done < "$LISTS/west_ds.txt"

# --- combined downsampled cohort (east + downsampled west) ---
cat "$LISTS/east.txt" "$LISTS/west_ds.txt" > "$LISTS/downsample.txt"
cat "$LISTS/east_bams.txt" "$LISTS/west_ds_bams.txt" > "$LISTS/downsample_bams.txt"

echo "Excluded samples:"
printf "  %s\n" "${!EXCLUDE[@]}"

echo "Full cohort: $(wc -l < "$PROJ/bams_list.txt")"
echo "East: $NE   West: $NW   West downsampled -> $(wc -l < "$LISTS/west_ds.txt")"
echo "Downsampled cohort: $(wc -l < "$LISTS/downsample.txt") (seed=$SEED)"

# sanity check: leading zeros intact
echo "Leading-zero check:"
grep -E '/0[0-9]+\.bam$' "$PROJ/bams_list.txt" | head -n 3 || true
