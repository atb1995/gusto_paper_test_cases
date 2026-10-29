#!/bin/bash
#SBATCH --job-name=siqn_weak_scaling
#SBATCH --output=logs/siqn_weak_scaling_%j.out
#SBATCH --time=01:00:00
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
RESULTS_CSV=/projects/u35z/gusto_paper_test_cases/results/siqn_weak_scaling_underpop_oct26.csv
PETSC_OPTS_FILE=${SCRIPT_DIR}/fieldsplit0_opts.txt
DT=${DT:-300}
TMAX=${TMAX:-86400}
NCELLS=${NCELLS}   # must be provided, e.g. sbatch --export=ALL,NCELLS=64
NLAYERS=15

# -----------------------------------------
if [[ -n "${SLURM_NTASKS:-}" ]]; then
NPROCS="${SLURM_NTASKS}"
else
NPROCS=$(( ${SLURM_NNODES:-1} * ${SLURM_NTASKS_PER_NODE:-${NTASKS_PER_NODE:-1}} ))
fi
DIRNAME="baro_weak_underpop_dt${DT}_C${NCELLS}_L${NLAYERS}_n${NPROCS}"
mkdir -p logs
mkdir -p "$(dirname "${RESULTS_CSV}")"

echo "=== [Starting weak scaling: C${NCELLS}, L${NLAYERS}, dt=${DT}, tmax=${TMAX}, nprocs=${NPROCS}] ==="
echo ">>> dirname=${DIRNAME}"
echo ">>> cache=${CACHE_ROOT}"
date

rm -rf "/projects/u35z/gusto_paper_test_cases/batch_scripts/isambard/siqn_scaling/results/${DIRNAME}"
rm -rf "${SCRIPT_DIR}/results/${DIRNAME}"

# -----------------------------------------
# Launch command for parallel runs
# -----------------------------------------
launch_cmd="srun --kill-on-bad-exit=1 --hint=nomultithread --cpu-bind=cores --unbuffered"

# -----------------------------------------
# Warm the kernel cache: tiny horizontal mesh, SAME number of layers
# -----------------------------------------
srun -N 1 -n 1 python3 "${PYTHON_SCRIPT}" \
--dt "${DT}" \
--tmax "${DT}" \
--ncell_per_edge 4 \
--nlayers "${NLAYERS}" \
--dirname "${DIRNAME}_warmup" \
-petscpartitioner_type ptscotch -petscpartitioner_view :part_view.log \
-options_file "${PETSC_OPTS_FILE}"
echo "=== [Cache warm-up complete] ==="
date

# -----------------------------------------
# Report process placement
# -----------------------------------------
${launch_cmd} ${HOME}/xthi/src/xthi_mpi_mp 2>&1 | tee "logs/xthi_C${NCELLS}_n${NPROCS}.log"
echo "=== [Placement check complete] ==="
date

# -----------------------------------------
# Run
# -----------------------------------------
${launch_cmd} python3 "${PYTHON_SCRIPT}" \
--dt "${DT}" \
--tmax "${TMAX}" \
--ncell_per_edge "${NCELLS}" \
--nlayers "${NLAYERS}" \
--dirname "${DIRNAME}" \
--results_csv "${RESULTS_CSV}" \
--io_on False \
-log_view :logview_C${NCELLS}_n${NPROCS}.txt \
-petscpartitioner_type ptscotch -petscpartitioner_view :part_view.log \
-options_file "${PETSC_OPTS_FILE}" \
2>&1 | tee "logs/siqn_weak_scaling_C${NCELLS}_n${NPROCS}.log"

echo "=== [Run complete] ==="
date