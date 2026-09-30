#!/bin/bash
# ==============================================================================
# AutoPas Functor Benchmark Parameter Sweep Runner
# ==============================================================================
# Can be run locally or within an interactive cluster allocation (e.g. salloc).
#
# Examples:
#   # Default sweep across ATM & ATM2 on all kernels:
#   ./cluster/run_sweep.sh
#
#   # Quick sanity sweep:
#   ./cluster/run_sweep.sh --quick
#
#   # Custom parameters:
#   ./cluster/run_sweep.sh -f ATM,ATM2 -k SoATriple -p 16,32,64 --n3 both
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Find executable across common build directories
POSSIBLE_DIRS=(
    "${BUILD_DIR:-}"
    "${REPO_DIR}/cmake-build-relwithdebinfo-wsl---gcc"
    "${REPO_DIR}/cmake-build-release"
    "${REPO_DIR}/cmake-build-relwithdebinfo"
    "${REPO_DIR}/build"
)
BIN=""
for dir in "${POSSIBLE_DIRS[@]}"; do
    if [ -n "${dir}" ] && [ -f "${dir}/AP_Functor_Bench" ]; then
        BIN="${dir}/AP_Functor_Bench"
        BUILD_DIR="${dir}"
        break
    fi
done

if [ -z "${BIN}" ]; then
    echo "Error: AP_Functor_Bench binary not found."
    echo "Please build the project first."
    exit 1
fi

# Defaults
FUNCTORS="ATM,ATM2"
KERNELS="SoASingle,SoAPair,SoATriple,AoS"
PARTICLES="16,32,64,128,256"
NEWTON3="both"
POOL_SIZE="1000"
REPETITIONS="1"
MIN_TIME="0.2s"
QUICK=false

# Parse command line overrides
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --quick)
            QUICK=true
            PARTICLES="16,32"
            POOL_SIZE="10"
            REPETITIONS="1"
            MIN_TIME="0.01s"
            shift
            ;;
        -f|--functor)
            FUNCTORS="$2"
            shift 2
            ;;
        -k|--kernel)
            KERNELS="$2"
            shift 2
            ;;
        -p|--particles)
            PARTICLES="$2"
            shift 2
            ;;
        --n3)
            NEWTON3="$2"
            shift 2
            ;;
        -l|--label|--tag)
            LABEL="$2"
            shift 2
            ;;
        --pool-size)
            POOL_SIZE="$2"
            shift 2
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

LABEL="${LABEL:-}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
if [ -n "${LABEL}" ]; then
    OUT_DIR="${REPO_DIR}/results/${TIMESTAMP}_${LABEL}"
else
    OUT_DIR="${REPO_DIR}/results/sweep_${TIMESTAMP}"
fi
mkdir -p "${OUT_DIR}"
JSON_FILE="${OUT_DIR}/results.json"

echo "=========================================="
echo "Starting Functor Benchmark Sweep"
echo "Binary:     ${BIN}"
echo "Functors:   ${FUNCTORS}"
echo "Kernels:    ${KERNELS}"
echo "Particles:  ${PARTICLES}"
echo "Newton-3:   ${NEWTON3}"
echo "Pool Size:  ${POOL_SIZE}"
echo "Reps:       ${REPETITIONS}"
echo "Min Time:   ${MIN_TIME}"
echo "Output:     ${OUT_DIR}"
echo "=========================================="

# OpenMP single-thread affinity
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

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
    --benchmark_out_format=json \
    ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}

if [ -n "${LABEL}" ]; then
    ln -sf "${JSON_FILE}" "${REPO_DIR}/results/${LABEL}.json"
    echo "Created convenient shortcut: results/${LABEL}.json -> ${JSON_FILE}"
fi

echo ""
echo "Benchmark completed. Running automated analysis..."

if command -v python3 &>/dev/null && [ -f "${REPO_DIR}/scripts/analyze_benchmarks.py" ]; then
    python3 "${REPO_DIR}/scripts/analyze_benchmarks.py" "${JSON_FILE}" \
        --baseline ATM \
        --output-dir "${OUT_DIR}/plots"
    echo "Plots and tables saved to: ${OUT_DIR}/plots"
fi

echo "=========================================="
echo "Sweep Complete: ${OUT_DIR}"
echo "=========================================="
