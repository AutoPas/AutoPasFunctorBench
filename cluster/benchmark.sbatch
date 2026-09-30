#!/bin/bash
#SBATCH --job-name=ap_functor_bench
#SBATCH --output=results/%j_benchmark.log
#SBATCH --error=results/%j_benchmark.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --exclusive

# Stop on errors
set -euo pipefail

# -------------------------------------------------------------
# Environment & Module Setup
# -------------------------------------------------------------
# module --force purge
# module load gcc/12 cmake/3.26 python/3.10

# Disable multi-threading / OpenMP to ensure single-core evaluation
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_DISPLAY_AFFINITY=TRUE

# -------------------------------------------------------------
# Paths & Directories
# -------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${REPO_DIR}/build}"
BIN="${BUILD_DIR}/AP_Functor_Bench"

JOB_ID="${SLURM_JOB_ID:-manual}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LABEL="${LABEL:-}"
if [ -n "${LABEL}" ]; then
    OUT_DIR="${REPO_DIR}/results/${TIMESTAMP}_job${JOB_ID}_${LABEL}"
else
    OUT_DIR="${REPO_DIR}/results/${TIMESTAMP}_job${JOB_ID}"
fi
mkdir -p "${OUT_DIR}"

JSON_FILE="${OUT_DIR}/benchmark_results.json"
NODE_INFO="${OUT_DIR}/node_info.txt"

echo "=========================================="
echo "AutoPas Functor Benchmark Cluster Run"
echo "Job ID:    ${JOB_ID}"
echo "Node:      $(hostname)"
echo "Timestamp: ${TIMESTAMP}"
echo "Output:    ${OUT_DIR}"
echo "=========================================="

# -------------------------------------------------------------
# Capture Hardware & System Metadata
# -------------------------------------------------------------
{
    echo "Hostname: $(hostname)"
    echo "Kernel:   $(uname -r)"
    echo "Date:     $(date)"
    echo "------------------ lscpu ------------------"
    lscpu 2>/dev/null || true
    echo "------------------ CPU Scaling Governor ------------------"
    cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || true
} > "${NODE_INFO}"

# -------------------------------------------------------------
# Check Executable
# -------------------------------------------------------------
if [ ! -f "${BIN}" ]; then
    echo "Error: Binary not found at ${BIN}."
    echo "Please build the project first:"
    echo "  cmake -B ${BUILD_DIR} -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=\"-march=native\""
    echo "  cmake --build ${BUILD_DIR} --target AP_Functor_Bench -j"
    exit 1
fi

# -------------------------------------------------------------
# Benchmark Configuration (Override via environment variables)
# -------------------------------------------------------------
FUNCTORS="${FUNCTORS:-ATM,ATM2}"
KERNELS="${KERNELS:-SoASingle,SoAPair,SoATriple,AoS}"
PARTICLES="${PARTICLES:-16,32,64,128,256}"
NEWTON3="${NEWTON3:-both}"
POOL_SIZE="${POOL_SIZE:-1000}"
REPETITIONS="${REPETITIONS:-1}"
MIN_TIME="${MIN_TIME:-0.5s}"

echo "Functors:   ${FUNCTORS}"
echo "Kernels:    ${KERNELS}"
echo "Particles:  ${PARTICLES}"
echo "Newton-3:   ${NEWTON3}"
echo "Pool Size:  ${POOL_SIZE}"
echo "Reps:       ${REPETITIONS}"
echo "Min Time:   ${MIN_TIME}"
echo "------------------------------------------"

# -------------------------------------------------------------
# Execute Benchmark
# -------------------------------------------------------------
"${BIN}" \
    -v \
    -f "${FUNCTORS}" \
    -k "${KERNELS}" \
    -p "${PARTICLES}" \
    --pool-size "${POOL_SIZE}" \
    --n3 "${NEWTON3}" \
    --benchmark_repetitions="${REPETITIONS}" \
    --benchmark_report_aggregates_only=true \
    --benchmark_min_time="${MIN_TIME}" \
    --benchmark_out="${JSON_FILE}" \
    --benchmark_out_format=json

echo "Benchmark finished successfully. Results saved to: ${JSON_FILE}"

if [ -n "${LABEL}" ]; then
    ln -sf "${JSON_FILE}" "${REPO_DIR}/results/${LABEL}.json"
    echo "Created convenient shortcut: results/${LABEL}.json -> ${JSON_FILE}"
fi

# -------------------------------------------------------------
# Post-Processing & Automated Analysis
# -------------------------------------------------------------
if command -v python3 &>/dev/null && [ -f "${REPO_DIR}/scripts/analyze_benchmarks.py" ]; then
    echo ""
    echo "Running analysis and generating plots..."
    python3 "${REPO_DIR}/scripts/analyze_benchmarks.py" "${JSON_FILE}" \
        --baseline ATM \
        --output-dir "${OUT_DIR}/plots"
    echo "Plots saved to: ${OUT_DIR}/plots"
fi

echo "=========================================="
echo "Job Complete: ${OUT_DIR}"
echo "=========================================="
