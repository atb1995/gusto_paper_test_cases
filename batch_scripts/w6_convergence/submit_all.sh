#!/bin/bash
# submit_w6.sh — queues one Slurm job per scheme per dt
# Usage: bash submit_w6.sh [explicit|implicit|all]
set -euo pipefail

SLURM_SCRIPT="$(dirname "$0")/run_w6.slurm"

IMPLICIT_SCHEMES=(siqn imex_ssp3 imex_ark4 backward_euler)
#IMPLICIT_SCHEMES=(backward_euler)
IMPLICIT_DTS=(480 240 120 60 30)

EXPLICIT_SCHEMES=(heun ssprk3)
EXPLICIT_DTS=(60 30 15 7.5 3.75 1.875)

MODE=${1:-all}

run_schemes() {
    local schemes=("${!1}")
    local dts=("${!2}")
    local total=$(( ${#schemes[@]} * ${#dts[@]} ))
    echo "Submitting ${total} jobs (${#schemes[@]} schemes x ${#dts[@]} dt values):"
    echo "  Schemes: ${schemes[*]}"
    echo "  DTs:     ${dts[*]}"
    echo ""
    for scheme in "${schemes[@]}"; do
        for dt in "${dts[@]}"; do
            job_id=$(sbatch \
                --export=ALL,SCHEME="${scheme}",DT="${dt}" \
                "${SLURM_SCRIPT}" | awk '{print $NF}')
            echo "  Submitted scheme=${scheme}, dt=${dt}  →  job ${job_id}"
        done
        echo ""
    done
}

case "$MODE" in
    explicit)
        run_schemes EXPLICIT_SCHEMES[@] EXPLICIT_DTS[@]
        ;;
    implicit)
        run_schemes IMPLICIT_SCHEMES[@] IMPLICIT_DTS[@]
        ;;
    all)
        run_schemes EXPLICIT_SCHEMES[@] EXPLICIT_DTS[@]
        run_schemes IMPLICIT_SCHEMES[@] IMPLICIT_DTS[@]
        ;;
    *)
        echo "ERROR: unknown mode '${MODE}'. Use explicit, implicit, or all." >&2
        exit 1
        ;;
esac

echo "Done. Check queue with: squeue -u \$USER"
