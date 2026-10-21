#!/usr/bin/env bash
# submit_rk_all.sh — queues one Slurm job per RK scheme
# Usage: bash submit_rk_all.sh

set -euo pipefail

SLURM_SCRIPT="$(dirname "$0")/run_rk_w6.slurm"

declare -A SCHEMES
#SCHEMES["EXP_RK"]="exp_rk"
SCHEMES["RK4"]="rk4"
#SCHEMES["RK5"]="rk5"

echo "Submitting ${#SCHEMES[@]} jobs:"
echo ""

for scheme in "${!SCHEMES[@]}"; do
    case_dir="${SCHEMES[$scheme]}"
    job_id=$(sbatch \
        --export=ALL,SCHEME="${scheme}",CASE_DIR="${case_dir}" \
        "${SLURM_SCRIPT}" | awk '{print $NF}')
    echo "  Submitted scheme=${scheme}, dir=test_cases/${case_dir}  →  job ${job_id}"
done

echo ""
echo "Done. Check queue with: squeue -u \$USER"

