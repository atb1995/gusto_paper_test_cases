#!/bin/bash
#SBATCH --job-name=siqn_dry_baro
#SBATCH --output=logs/siqn_dry_baro_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=72
#SBATCH --mem=220G
#SBATCH --distribution=block:block
#SBATCH --exclusive
set -euo pipefail
# -----------------------------------------
# Modules & environment
# -----------------------------------------
module load PrgEnv-gnu
module load cray-python
module swap cray-mpich/9.1.0 cray-mpich/8.1.30
export GUSTO_DIR=~/gusto_apr26
export GUSTO_ENV=/home/u35z/alexbrown.u35z/gusto_apr26/venv-gusto
export PETSC_DIR=~/gusto_apr26/petsc PETSC_ARCH=arch-firedrake-default
export LIBSUPERMESH_DIR=$GUSTO_ENV/lib/python3.11/site-packages/libsupermesh
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LIBSUPERMESH_DIR/lib:$LD_LIBRARY_PATH
export PYOP2_CACHE_DIR=/projects/u35z/pyop2_cache
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR=/projects/u35z/tsfc_cache
mkdir -p $PYOP2_CACHE_DIR $FIREDRAKE_TSFC_KERNEL_CACHE_DIR
export XDG_CACHE_HOME=/projects/u35z/cache
mkdir -p $XDG_CACHE_HOME
source $GUSTO_ENV/bin/activate
export PYTHONPATH=/home/u35z/alexbrown.u35z/gusto_apr26/dc_improvements:$PYTHONPATH
export OMP_NUM_THREADS=1
# -----------------------------------------
# User configuration
# -----------------------------------------
PYTHON_SCRIPT=/projects/u35z/gusto_paper_test_cases/test_cases/siqn/compressible_euler/dry_baroclinic_sphere.py
DT=${DT:-450}
TMAX=${TMAX:-864000}
NCELLS=96
NLAYERS=30
# -----------------------------------------
DIRNAME="dry_baro_dt${DT}_C${NCELLS}_L${NLAYERS}"
mkdir -p logs
rm -rf "/projects/u35z/gusto_paper_test_cases/results/${DIRNAME}"
echo "=== [Starting dry baroclinic sphere: C${NCELLS}, L${NLAYERS}, dt=${DT}, tmax=${TMAX}] ==="
echo ">>> dirname=${DIRNAME}"
date
# -----------------------------------------
# Run
# -----------------------------------------
srun --hint=nomultithread \
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

