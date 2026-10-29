#!/bin/bash
#SBATCH --job-name=siqn_dry_baro
#SBATCH --output=logs/siqn_dry_baro_%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=72
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
PYTHON_SCRIPT=/projects/u35z/gusto_paper_test_cases/test_cases/siqn/compressible_euler/dry_baroclinic_sphere.py
DT=${DT:-450}
TMAX=${TMAX:-864000}
NCELLS=96
NLAYERS=15
# -----------------------------------------
DIRNAME="dry_baro_dt${DT}_C${NCELLS}_L${NLAYERS}"
mkdir -p logs
rm -rf "/projects/u35z/gusto_paper_test_cases/results/${DIRNAME}"
echo "=== [Starting dry baroclinic sphere: C${NCELLS}, L${NLAYERS}, dt=${DT}, tmax=${TMAX}] ==="
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
--nlayers 10 \
--dirname "${DIRNAME}_warmup"
echo "=== [Cache warm-up complete] ==="
date

# -----------------------------------------
# Run
# -----------------------------------------
srun --kill-on-bad-exit=1 \
--hint=nomultithread \
--cpu_bind=cores \
python3 "${PYTHON_SCRIPT}" \
--dt "${DT}" \
--tmax "${TMAX}" \
--ncell_per_edge "${NCELLS}" \
--nlayers "${NLAYERS}" \
--dumpfreq 192 \
--dirname "${DIRNAME}" \
2>&1 | tee "logs/siqn_dry_baro_dt${DT}.log"
echo "=== [Run complete] ==="
date
