#!/bin/bash
# Submit SIQN strong scaling jobs: 8, 16, 32, 64, 128, 256 cores

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/dryrun_strong_scaling.slurm"

# ntasks_per_node:nodes:walltime
declare -A CORE_CONFIG
CORE_CONFIG[8]="8:1:08:00:00"
CORE_CONFIG[16]="16:1:04:00:00"
CORE_CONFIG[32]="32:1:02:00:00"
CORE_CONFIG[64]="64:1:01:30:00"
CORE_CONFIG[128]="128:1:01:00:00"
CORE_CONFIG[256]="128:2:00:30:00"

for NCORES in 8 16 32 64 128 256; do
    IFS=':' read -r NTASKS NNODES WALLTIME <<< "${CORE_CONFIG[$NCORES]}"
    echo "Submitting: ${NCORES} cores (${NNODES} node(s), ${NTASKS} tasks/node, time=${WALLTIME})"
    sbatch \
        --nodes="${NNODES}" \
        --ntasks-per-node="${NTASKS}" \
        --time="${WALLTIME}" \
        --export=ALL,NTASKS_PER_NODE="${NTASKS}" \
        "${SLURM_SCRIPT}"
done
