#!/bin/bash
#SBATCH --job-name=siqn_dryrun_strong
#SBATCH --output=logs/siqn_dryrun_strong_%j.out
#SBATCH --mem=220G
#SBATCH --distribution=block:block
#SBATCH --exclusive
set -euo pipefail

# -----------------------------------------
# Modules & environment
# -----------------------------------------
module load PrgEnv-gnu
module load cray-python
export GUSTO_DIR=~/gusto_oct26
export GUSTO_ENV=/home/u35z/alexbrown.u35z/gusto_oct26/venv-gusto
export PETSC_DIR=/home/u35z/alexbrown.u35z/gusto_oct26/petsc
export PETSC_ARCH=arch-firedrake-default
export LIBSUPERMESH_DIR=$GUSTO_ENV/lib/python3.11/site-packages/libsupermesh
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LIBSUPERMESH_DIR/lib:$LD_LIBRARY_PATH

# Per-job caches on shared storage (visible to all nodes, isolated from other jobs)
CACHE_ROOT=/projects/u35z/gusto_paper_test_cases/job_caches/${SLURM_JOB_ID}
export PYOP2_CACHE_DIR=$CACHE_ROOT/pyop2
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR=$CACHE_ROOT/tsfc
export XDG_CACHE_HOME=$CACHE_ROOT/xdg
mkdir -p "$PYOP2_CACHE_DIR" "$FIREDRAKE_TSFC_KERNEL_CACHE_DIR" "$XDG_CACHE_HOME"

export FI_CXI_RX_MATCH_MODE=software

source $GUSTO_ENV/bin/activate
export OMP_NUM_THREADS=1
# -----------------------------------------
# User configuration
# -----------------------------------------
SCRIPT_DIR=/projects/u35z/gusto_paper_test_cases/test_cases/siqn/compressible_euler
PYTHON_SCRIPT=${SCRIPT_DIR}/dry_baroclinic_sphere.py
RESULTS_CSV=/projects/u35z/gusto_paper_test_cases/results/siqn_dryrun_strong_scaling_oct26.csv
DT=${DT:-225}
TMAX=${TMAX:-21600}
NCELLS=96
NLAYERS=15

# -----------------------------------------
# Prefer Slurm's computed global task count; fallback supports manual dry-runs.
if [[ -n "${SLURM_NTASKS:-}" ]]; then
NPROCS="${SLURM_NTASKS}"
else
NPROCS=$(( ${SLURM_NNODES:-1} * ${SLURM_NTASKS_PER_NODE:-${NTASKS_PER_NODE:-1}} ))
fi
DIRNAME="dry_baro_strong_dt${DT}_C${NCELLS}_L${NLAYERS}_n${NPROCS}"
mkdir -p logs
mkdir -p "$(dirname "${RESULTS_CSV}")"

echo "=== [Starting dry baroclinic sphere: C${NCELLS}, L${NLAYERS}, dt=${DT}, tmax=${TMAX}, nprocs=${NPROCS}] ==="
echo ">>> dirname=${DIRNAME}"
echo ">>> cache=${CACHE_ROOT}"
date

# -----------------------------------------
# Warm the kernel cache: tiny horizontal mesh, SAME number of layers
# -----------------------------------------
srun -N 1 -n 1 python3 "${PYTHON_SCRIPT}" \
--dt "${DT}" \
--tmax "${DT}" \
--ncell_per_edge 4 \
--nlayers "${NLAYERS}" \
--dirname "${DIRNAME}_warmup"
echo "=== [Cache warm-up complete] ==="
date

# -----------------------------------------
# Run
# -----------------------------------------
rm -rf "/projects/u35z/gusto_paper_test_cases/results/${DIRNAME}"
rm -rf "${SCRIPT_DIR}/results/${DIRNAME}"
srun --kill-on-bad-exit=1 \
--hint=nomultithread \
--cpu-bind=cores \
python3 "${PYTHON_SCRIPT}" \
--dt "${DT}" \
--tmax "${TMAX}" \
--ncell_per_edge "${NCELLS}" \
--nlayers "${NLAYERS}" \
--dirname "${DIRNAME}" \
--results_csv "${RESULTS_CSV}" \
--io_on False \
2>&1 | tee "logs/siqn_dryrun_strong_dt${DT}_n${NPROCS}.log"

echo "=== [Run complete] ==="
date