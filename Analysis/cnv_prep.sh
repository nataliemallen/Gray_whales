#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=cnv_prep
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# one-time cnvpytor setup: gc-content file + gray_whale reference config (run before cnvs.sh)


set -euo pipefail
module --force purge
module load biocontainers
module load samtools
module load anaconda
conda activate cnvpytor
 
PROJ="/scratch/gautschi/allen715/2026_whales"
REF="$PROJ/reference/ref.fa"
FAI="$REF.fai"
CHRS="$PROJ/reference/chrs.txt"
OUT="$PROJ/cnvs"
GCPYTOR="$OUT/gray_whale_gc_file.pytor"
CONF="$OUT/whale_config.py"
 
mkdir -p "$OUT/pytor_files" "$OUT/results"
cd "$OUT"
[ -f "$FAI" ] || samtools faidx "$REF"
 
# 1) GC-content file (keyed by the reference's chromosome names)
echo "[$(date)] building GC file ..."
cnvpytor -root "$GCPYTOR" -gc "$REF" -make_gc_file
 
# 2) auto-generate the reference-genome config from the .fai + chrs.txt
echo "[$(date)] writing $CONF ..."
python3 - "$FAI" "$CHRS" "$GCPYTOR" "$CONF" <<'PY'
import sys
from collections import OrderedDict
fai, chrs_file, gc_file, out = sys.argv[1:5]
 
lengths = {}
for line in open(fai):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 2:
        lengths[f[0]] = int(f[1])          # name -> length (exact string, keeps zeros)
 
chroms = OrderedDict()
missing = []
for line in open(chrs_file):
    c = line.strip()
    if not c:
        continue
    if c in lengths:
        chroms[c] = (lengths[c], "A")      # "A" = autosome-type; fine for all here
    else:
        missing.append(c)
 
# CNVpytor loads the config by exec() and then iterates a dict named
# `import_reference_genomes` (genome.py: `for g in import_reference_genomes`).
# The config MUST define that variable -- assigning to
# cnvpytor.Genome.reference_genomes[...] does NOT work.
with open(out, "w") as fh:
    fh.write("from collections import OrderedDict\n\n")
    fh.write("import_reference_genomes = {\n")
    fh.write('    "gray_whale": {\n')
    fh.write('        "name": "gray_whale",\n')
    fh.write('        "species": "Eschrichtius robustus",\n')
    fh.write('        "chromosomes": OrderedDict([\n')
    for c, (l, t) in chroms.items():
        fh.write(f'            ({c!r}, ({l}, {t!r})),\n')
    fh.write("        ]),\n")
    fh.write(f'        "gc_file": {gc_file!r},\n')
    fh.write("    }\n")
    fh.write("}\n")
 
print(f"  chromosomes registered: {len(chroms)}")
if missing:
    print(f"  WARNING: {len(missing)} chrs.txt names not found in .fai (skipped): {missing[:5]}")
PY
 
echo "[$(date)] done. Config: $CONF"
echo "Next: debug one sample with   bash cnv.sh <SAMPLE_ID>   then submit the array."

