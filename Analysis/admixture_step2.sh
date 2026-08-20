#!/bin/bash
#SBATCH -A dewoody
#SBATCH --job-name=admix_eval
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 96
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-user=allen715@purdue.edu
#SBATCH --mail-type=END,FAIL
# admixture model selection: log-likelihoods, deltaK, evaladmix (Table S1, Fig S7)

module --force purge
module load biocontainers
module load angsd

PROJ="/scratch/gautschi/allen715/2026_whales"
ADX="$PROJ/admixture/ADX"
beagle="$PROJ/PCA/whale.beagle.gz"
cd "$PROJ/admixture"

# # 1) log-likelihoods: "<run> <loglike>"
# > log_likelihoods.txt
# for log_file in "$ADX"/*.log; do
#   run_name=$(grep "outfiles=" "$log_file" | sed 's/.*outfiles=//')
#   log_like=$(grep "best like=" "$log_file" | sed 's/.*best like=//;s/ after.*//')
#   echo "$run_name $log_like" >> log_likelihoods.txt
# done
# 
# # 2) Evanno deltaK (writes deltaK_results.csv)
# module load anaconda 2>/dev/null || true
# python3 "$PROJ/calc_deltaK.py"

# 3) evalAdmix on the best (highest-likelihood) run per K
EVAL=/home/allen715/evalAdmix/evalAdmix
for K in 1 2 3 4 5; do
  best=$(grep "admix_K${K}_run" log_likelihoods.txt | sort -k2,2g | tail -n1 | awk '{print $1}')
  echo "K=$K best run: $best"
  $EVAL -beagle "$beagle" -fname "$ADX/${best}.fopt.gz" -qname "$ADX/${best}.qopt" -P 32
  mv output.corres.txt "K${K}_output.corres.txt"
  cp "$ADX/${best}.qopt" "admix_K${K}_best.qopt"   # for eval_admix.R / admix_plot.R
done
