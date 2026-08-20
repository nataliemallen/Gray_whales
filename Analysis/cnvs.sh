#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=cnv
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%A_%a.err
#SBATCH -o %x_%A_%a.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-64%16
# cnv calling with cnvpytor at 10 kb and 100 kb, downsampled dataset (Fig 5A,B, S9)
# requires cnv_prep.sh; run: sbatch cnvs.sh (array) or bash cnvs.sh <sample_id> (debug)

module --force purge
module load biocontainers
module load samtools
module load anaconda

# make conda activate work in a non-interactive sbatch shell
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cnvpytor

# fail fast if the env didn't activate (prevents silent empty output)
if ! command -v cnvpytor >/dev/null 2>&1; then
  echo "ERROR: cnvpytor not found after 'conda activate cnvpytor'." >&2
  echo "  python: $(command -v python)   conda base: $(conda info --base 2>/dev/null)" >&2
  exit 1
fi
echo "using cnvpytor: $(command -v cnvpytor)"

PROJ="/scratch/gautschi/allen715/2026_whales"
BAMDIR="$PROJ/bams"
OUT="$PROJ/cnvs"
CONF="$OUT/whale_config.py"
RG="gray_whale"
USE_GC=${USE_GC:-1}     # 1 = GC-corrected calling via the custom gray whale genome
                        #     (recommended; requires cnv_prep.sh to have been run)
                        # 0 = no custom genome / no GC correction (fallback)
cd "$OUT"

# config sanity (gc-corrected mode only)
if [ "$USE_GC" = 1 ] && [ ! -f "$CONF" ]; then
  echo "ERROR: $CONF not found. Run cnv_prep.sh first, or set USE_GC=0." >&2
  exit 1
fi

# --- pick sample + mode ---
if [ -n "${1:-}" ]; then
  ID="$1"; DEBUG=1
  echo "### DEBUG MODE: single sample $ID ###"
elif [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJ/lists/downsample.txt"); DEBUG=0
else
  echo "Usage: sbatch cnv.sh   (production array)   |   bash cnv.sh <SAMPLE_ID>   (debug)" >&2
  exit 1
fi
[ -z "$ID" ] && { echo "No sample ID resolved." >&2; exit 1; }

BAM="${BAMDIR}/${ID}.bam"
PYTOR="${OUT}/pytor_files/${ID}.pytor"
[ -f "$BAM" ]      || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[ -f "${BAM}.bai" ] || samtools index "$BAM"

echo "[$(date)] CNVpytor start: $ID"

# optional debug: check bam @SQ names match the config chromosomes
if [ "$DEBUG" = 1 ] && [ "$USE_GC" = 1 ]; then
  echo "--- BAM @SQ names (first 5) ---"
  samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); print $2}' | head -5
  echo "--- config chromosomes (first 5) ---"
  grep -oE "\('[^']+'" "$CONF" | tr -d "('" | head -5
  echo "(these should overlap; if not, calls will be empty)"
fi

# start from a clean root (a stale .pytor can lack gc/reference and misbehave)
rm -f "$PYTOR"

# -conf flag is only used in GC-corrected mode (path has no spaces, so unquoted)
CONF_FLAG=""; [ "$USE_GC" = 1 ] && CONF_FLAG="-conf $CONF"

# 1) import read depth
cnvpytor $CONF_FLAG -root "$PYTOR" -rd "$BAM"

# guard: if -rd was killed (e.g. oom) the root is missing/empty
if [ ! -s "$PYTOR" ]; then
  echo "ERROR: read-depth import failed (${PYTOR} missing/empty)." >&2
  echo "       Most common cause is an OOM kill -- raise --mem in the SBATCH header." >&2
  exit 1
fi

# 2) register genome (gc mode only) -> histograms -> partition -> call
if [ "$USE_GC" = 1 ]; then
  cnvpytor -conf "$CONF" -root "$PYTOR" -rg "$RG"
fi
cnvpytor $CONF_FLAG -root "$PYTOR" -his 10000 100000
cnvpytor $CONF_FLAG -root "$PYTOR" -partition 10000 100000
cnvpytor $CONF_FLAG -root "$PYTOR" -call 10000  > "${OUT}/results/${ID}.10kb.calls.tsv"
cnvpytor $CONF_FLAG -root "$PYTOR" -call 100000 > "${OUT}/results/${ID}.100kb.calls.tsv"

# --- report / debug diagnostics ---
if [ "$DEBUG" = 1 ]; then
  echo "--- signals stored in ${ID}.pytor ---"
  cnvpytor -root "$PYTOR" -ls || true
fi
for b in 10kb 100kb; do
  f="${OUT}/results/${ID}.${b}.calls.tsv"
  n=$(wc -l < "$f")
  echo "  $f : $n calls"
  if [ "$n" -eq 0 ]; then
    echo "  WARNING: 0 calls at ${b}. Check that cnv_prep.sh ran, that BAM @SQ names"
    echo "           match the config chromosomes, and that read depth is sufficient."
  fi
done

echo "[$(date)] CNVpytor done: $ID"