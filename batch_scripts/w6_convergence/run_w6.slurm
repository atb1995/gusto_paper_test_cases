#!/bin/bash
#SBATCH --job-name=w6_convergence
#SBATCH --output=logs/w6_%j.out
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=60G
#SBATCH --exclusive
#SBATCH --distribution=block:block
set -euo pipefail
# -----------------------------------------
# Modules & environment
# -----------------------------------------
module purge
module load PrgEnv-gnu
module load cray-python
export GUSTO_ENV=/home/u35z/alexbrown.u35z/gusto_oct26/venv-gusto
export PYTHONPATH=/home/u35z/alexbrown.u35z/gusto_oct26/dc_improvements:$PYTHONPATH
export PETSC_DIR=/home/u35z/alexbrown.u35z/gusto_oct26/petsc PETSC_ARCH=arch-firedrake-default
export LIBSUPERMESH_DIR=$GUSTO_ENV/lib/python3.11/site-packages/libsupermesh
export LD_LIBRARY_PATH=/opt/cray/libfabric/1.22.0/lib64:$PETSC_DIR/$PETSC_ARCH/lib:$LIBSUPERMESH_DIR/lib:${LD_LIBRARY_PATH:-}
export PYOP2_CACHE_DIR=/projects/u35z/gusto_paper_test_cases/pyop2_cache
export XDG_CACHE_HOME=/projects/u35z/gusto_paper_test_cases/cache
source $GUSTO_ENV/bin/activate
export OMP_NUM_THREADS=1
# -----------------------------------------
# User configuration
# -----------------------------------------
NCELLS=128
TMAX=86400
CORES=24
# -----------------------------------------
# Validate injected variables
# -----------------------------------------
for var in SCHEME DT; do
    if [ -z "${!var:-}" ]; then
        echo "ERROR: ${var} not set. Use submit_w6.sh to submit this script." >&2
        exit 1
    fi
done
PYTHON_SCRIPT="/projects/u35z/gusto_paper_test_cases/test_cases/${SCHEME}/shallow_water/williamson_6.py"
dt="${DT}"
dumpfreq=$(python3 -c "print(int($TMAX / $dt))")
dirname="williamson6_dt$(echo ${dt} | sed 's/\./p/')_C${NCELLS}_${SCHEME}_j${SLURM_JOB_ID}"
# -----------------------------------------
# Prepare output directory
# -----------------------------------------
OUTPUT_DIR="/projects/u35z/gusto_paper_test_cases/test_cases/${SCHEME}/shallow_water"
mkdir -p "${OUTPUT_DIR}/logs"
cd "${OUTPUT_DIR}"
rm -rf "results/${dirname}"*
RESULTS_CSV="results_time_vs_error_${SCHEME}_C${NCELLS}.csv"
if [ ! -f "${RESULTS_CSV}" ]; then
    echo "scheme,dt_s,tmax_s,ncells,elapsed_s" > "${RESULTS_CSV}"
fi
OUTLOG="logs/${SCHEME}_C${NCELLS}_dt${dt}_j${SLURM_JOB_ID}.log"
echo "=== [Starting W6 convergence: scheme=${SCHEME}, C${NCELLS}, dt=${dt}, cores=${CORES}] ==="
date
echo ">>> dirname=${dirname}, dumpfreq=${dumpfreq}"
# -----------------------------------------
# Run
# -----------------------------------------
SCHEME="${SCHEME}" \
NCELLS="${NCELLS}" \
DT="${dt}" \
TMAX="${TMAX}" \
DIRNAME="${dirname}" \
srun --ntasks="${CORES}" \
    --hint=nomultithread \
    --cpu_bind=cores \
    python3 "${PYTHON_SCRIPT}" \
        --ncells_per_edge "${NCELLS}" \
        --dt "${dt}" \
        --tmax "${TMAX}" \
        --dumpfreq "${dumpfreq}" \
        --dirname "${dirname}" \
    2>&1 | tee "${OUTLOG}"
# Extract the RESULT line and append without the "RESULT," prefix
if grep -q "^RESULT," "${OUTLOG}"; then
    grep "^RESULT," "${OUTLOG}" | sed 's/^RESULT,//' >> "${RESULTS_CSV}"
else
    echo "WARNING: no RESULT line found in ${OUTLOG}" >&2
fi
echo "=== [Run complete] ==="
date
