#!/bin/bash
# Submit SIQN strong scaling jobs: 8, 16, 32, 64, 128, 256 cores

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/strong_scaling_paper.slurm"

# ntasks_per_node:nodes:walltime
declare -A NODE_CONFIG
# NODE_CONFIG[1]="1:01:00:00"
# NODE_CONFIG[2]="2:01:00:00"
NODE_CONFIG[4]="4:01:00:00"
NODE_CONFIG[8]="8:00:45:00"
NODE_CONFIG[16]="16:00:30:00"
NODE_CONFIG[32]="32:00:30:00"
NODE_CONFIG[64]="64:00:30:00"
NTASKS_PER_NODE=36
for NODES in 4 8 16 32 64; do
    IFS=':' read -r NNODES WALLTIME <<< "${NODE_CONFIG[$NODES]}"
    echo "Submitting: ${NNODES} node(s), time=${WALLTIME}"
    sbatch \
        --ntasks-per-node="${NTASKS_PER_NODE}" \
        --nodes="${NNODES}" \
        --time="${WALLTIME}" \
        "${SLURM_SCRIPT}"
done
