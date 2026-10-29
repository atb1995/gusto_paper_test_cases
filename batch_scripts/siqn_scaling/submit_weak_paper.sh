#!/bin/bash
# Submit weak scaling jobs: C8/2, C16/8, C32/32, C64/128, C128/512

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/weak_scaling_paper.slurm"

# ncells : ntasks_per_node : nodes
declare -A WEAK_CONFIG
WEAK_CONFIG[24]="36:1"
WEAK_CONFIG[48]="36:4"
WEAK_CONFIG[96]="36:16"
WEAK_CONFIG[192]="36:64"
#WEAK_CONFIG[384]="144:64"

# Per-resolution time step and run length.
declare -A DT_BY_NCELLS
DT_BY_NCELLS[24]="900"
DT_BY_NCELLS[48]="450"
DT_BY_NCELLS[96]="225"
DT_BY_NCELLS[192]="112.5"
DT_BY_NCELLS[384]="56.25"

declare -A TMAX_BY_NCELLS
TMAX_BY_NCELLS[24]="86400"
TMAX_BY_NCELLS[48]="43200"
TMAX_BY_NCELLS[96]="21600"
TMAX_BY_NCELLS[192]="10800"
TMAX_BY_NCELLS[384]="5400"

for NCELLS in $(printf '%s\n' "${!WEAK_CONFIG[@]}" | sort -n); do
    CONFIG="${WEAK_CONFIG[$NCELLS]:-}"
    if [[ -z "${CONFIG}" ]]; then
        echo "Error: no WEAK_CONFIG entry for NCELLS=${NCELLS}" >&2
        exit 1
    fi

    IFS=':' read -r NTASKS NNODES <<< "${CONFIG}"
    DT="${DT_BY_NCELLS[$NCELLS]:-300}"
    TMAX="${TMAX_BY_NCELLS[$NCELLS]:-86400}"
    NPROCS=$((NNODES * NTASKS))
    echo "Submitting: C${NCELLS}, ${NPROCS} cores (${NNODES} node(s), ${NTASKS} tasks/node), dt=${DT}, tmax=${TMAX}"
    sbatch \
        --nodes="${NNODES}" \
        --ntasks-per-node="${NTASKS}" \
        --export=ALL,NTASKS_PER_NODE="${NTASKS}",NCELLS="${NCELLS}",DT="${DT}",TMAX="${TMAX}" \
        "${SLURM_SCRIPT}"
done
