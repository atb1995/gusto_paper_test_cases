#!/bin/bash
# Submit weak scaling jobs: C8/2, C16/8, C32/32, C64/128, C128/512

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/weak_scaling.slurm"

# ncells : ntasks_per_node : nodes
declare -A WEAK_CONFIG
WEAK_CONFIG[8]="2:1"
WEAK_CONFIG[16]="8:1"
WEAK_CONFIG[32]="32:1"
WEAK_CONFIG[64]="64:2"
WEAK_CONFIG[128]="64:8"

for NCELLS in 8 16 32 64 128; do
    IFS=':' read -r NTASKS NNODES <<< "${WEAK_CONFIG[$NCELLS]}"
    NPROCS=$((NNODES * NTASKS))
    echo "Submitting: C${NCELLS}, ${NPROCS} cores (${NNODES} node(s), ${NTASKS} tasks/node)"
    sbatch \
        --nodes="${NNODES}" \
        --ntasks-per-node="${NTASKS}" \
        --export=ALL,NTASKS_PER_NODE="${NTASKS}",NCELLS="${NCELLS}" \
        "${SLURM_SCRIPT}"
done
