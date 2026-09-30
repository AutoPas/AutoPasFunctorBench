# AutoPas Functor Bench

A standalone benchmarking and verification harness for evaluating and optimizing AutoPas interaction functors (such as the 3-body Axilrod-Teller-Muto potential).

It allows side-by-side performance benchmarking and numerical correctness verification between a **Baseline** functor and one or more **Candidate** functor implementations across data layouts (AoS, SoA) and cell stencils.

---

## Features

- **Side-by-Side Functor Comparison**: Benchmark baseline and candidate functors in the same binary under identical cache and hardware conditions.
- **In-Memory Correctness Verification**: Check numerical agreement between baseline and candidate functors before or without benchmarking (`--verify`, `--verify-only`).
- **Flexible Execution Kernels**: Target specific kernels (`AoS`, `SoASingle`, `SoAPair`, `SoATriple`, or `all`).
- **Newton-3 Comparison in One Run**: Compare Newton-3 enabled vs. disabled side-by-side (`--n3 both`, `--n3 on`, `--n3 off`).
- **High-Resolution Throughput Metrics**: Reports interaction rate (`Interactions/s`), triplet evaluation rate (`Triplets/s`), time per triplet, and hit rates using Google Benchmark counters.
- **Configurable Benchmarking Ergonomics**: Support for explicit particle counts (`-p 8,16,32`), cell pool sizing (`--pool-size`), and deterministic seeding (`--seed`).

---

## Requirements

- **CMake**: >= 3.14
- **C++ Compiler**: Supporting C++20 (GCC 11+, Clang 13+)
- **Git**: Required by CMake's `FetchContent` to download AutoPas, Google Benchmark, and CLI11

---

## Build Instructions

```bash
# Configure the build
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build the benchmark executable
cmake --build build --target AP_Functor_Bench -j$(nproc)
```

To build with native vectorization flags on your machine:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native"
cmake --build build --target AP_Functor_Bench -j$(nproc)
```

---

## Usage

```bash
./build/AP_Functor_Bench [OPTIONS]
```

### CLI Options

| Flag | Description | Default |
|---|---|---|
| `-f, --functor` | Comma-separated list of functors to test (`ATM`, `ATM2`, `ATMGlobals`, `all`) | `all` |
| `-k, --kernel` | Comma-separated list of kernels (`AoS`, `SoASingle`, `SoAPair`, `SoATriple`, `all`) | `all` |
| `-p, --particles` | Explicit comma-separated list of particle counts (e.g. `-p 16,32,64`) | *(none)* |
| `--min`, `--max` | Particle range per cell when `-p` is not set (powers of 2) | `1`, `512` |
| `--n3` | Newton-3 configuration: `on`, `off`, or `both` | `on` |
| `--no-n3` | Shortcut to disable Newton-3 (`--n3 off`) | |
| `--pool-size` | Number of cell stencils pre-generated in memory pool | `1000` |
| `-c, --cell-size` | Simulation cell size | `3.0` |
| `-r, --cutoff` | Interaction cutoff radius | `3.0` |
| `-s, --seed` | Random seed for particle generation | `42` |
| `-v, --verify` | Verify numerical correctness against baseline before benchmarking | `false` |
| `--verify-only` | Run correctness verification and exit immediately | `false` |
| `--verify-baseline` | Functor to treat as baseline in verification | first from `-f` or `ATM` |
| `--verify-candidate` | Functor to treat as candidate in verification | second from `-f` or `ATM2` |
| `--verify-particles` | Number of particles per cell used in verification | `16` |
| `--verify-tol` | Maximum absolute force tolerance for verification pass | `1e-10` |

> Any unrecognized options (such as `--benchmark_filter`, `--benchmark_min_time`, `--benchmark_out`, etc.) are passed directly to Google Benchmark.

---

## Examples

### 1. Correctness Verification Only
Verify that candidate `ATM2` matches baseline `ATM` across all kernels:
```bash
./build/AP_Functor_Bench --verify-only
```
Or verify a specific kernel:
```bash
./build/AP_Functor_Bench --verify-only -k SoATriple
```

### 2. Side-by-Side Functor Comparison
Verify and compare `ATM` and `ATM2` side-by-side with Newton-3 both ON and OFF:
```bash
./build/AP_Functor_Bench -v -f ATM,ATM2 -k SoATriple -p 16,32 --n3 both
```

### 3. Fast Local Development Run
Use a small cell pool and short benchmark duration for instant iteration during kernel optimization:
```bash
./build/AP_Functor_Bench -k SoATriple -f ATM2 -p 16 --pool-size 10 --benchmark_min_time=0.01s
```

---

## Post-Processing & Plotting

Analyze Google Benchmark JSON outputs and generate publication-ready comparison plots:

```bash
# Analyze a benchmark JSON and plot candidate speedup relative to Baseline (ATM)
python3 scripts/analyze_benchmarks.py results.json --baseline ATM --output-dir plots/

# Or compare across multiple separate runs / branches:
python3 scripts/analyze_benchmarks.py \
  --baseline-file Baseline=results_base.json \
  --candidate "Unroll=results_iter1.json" \
  --candidate "AVX512=results_iter2.json" \
  --output-dir plots/
```

This generates:
- `speedup_vs_baseline.png`: Speedup bar chart relative to the $1.0\times$ baseline reference line.
- `throughput_scaling.png`: GigaTriplets/second scaling curves across particle counts.
- `time_per_triplet.png`: Time per triplet in picoseconds (lower is better).

An interactive Jupyter Notebook is also available at `plot_benchmark.ipynb`.

---

## Cluster Usage & Automated Sweeps

For running benchmarks on compute clusters or conducting automated parameter sweeps:

### 1. Interactive Sweeps (`cluster/run_sweep.sh`)
Run a sweep across functors, kernels, and particle counts with strict OpenMP core affinity and automatic post-processing:
```bash
# Run default sweep across ATM & ATM2
./cluster/run_sweep.sh

# Run quick sanity sweep
./cluster/run_sweep.sh --quick

# Custom sweep
./cluster/run_sweep.sh -f ATM,ATM2 -k SoATriple,SoAPair -p 16,32,64,128 --n3 both
```
Results and plots are automatically saved to `results/sweep_<timestamp>/`.

### 2. SLURM Batch Submission (`cluster/benchmark.sbatch`)
To submit a production benchmark job on a SLURM-managed cluster:
```bash
sbatch cluster/benchmark.sbatch
```
This script:
- Enforces single-thread execution (`OMP_NUM_THREADS=1`) and exclusive node usage.
- Records node hardware metadata (`lscpu`, CPU frequency governor, hostname).
- Configures standard OpenMP thread affinity (`OMP_PLACES=cores`, `OMP_PROC_BIND=spread`) respecting Slurm's CPU allocation.
- Saves results to `results/<timestamp>_job<jobid>/`.
- Automatically runs `scripts/analyze_benchmarks.py` upon completion to produce summary tables and plots.

---

## Vectorization & Roofline Analysis (Intel Advisor)

To deeply inspect vectorization efficiency, instruction mix, and Roofline model placement without noise from verification or container setup:

### 1. Build with Intel ITT API enabled
```bash
cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_ITT=ON
cmake --build build --target AP_Functor_Bench -j
```
*(When `-DENABLE_ITT=OFF` (the default), the codebase compiles with zero dependencies and zero overhead.)*

### 2. Run Automated Survey & Roofline Analysis
```bash
# Source Intel oneAPI / Advisor environment
source /opt/intel/oneapi/setvars.sh  # or: module load intel/advisor

# Run automated Survey + Trip Counts & FLOPs + HTML Roofline export
./scripts/run_advisor.sh -f ATM2 -k SoATriple -p 64
```
This produces:
- `advisor_results/roofline.html`: Standalone, interactive HTML Roofline model chart.
- Full Advisor database viewable in GUI: `advisor-gui advisor_results`.

---

## Developing New Functor Variants

To benchmark an optimized variant of an AutoPas functor:

1. **Add Functor Header**: Place your functor implementation in `functor-variants/` (e.g., `functor-variants/MyFunctor.h`).
2. **Register Functor**: In `main.cpp`, add registration in `initRegistry()`:
   ```cpp
   registry.registerFunctor<MyFunctor>(
       "MyFunctor",
       "Candidate Functor with SIMD optimization",
       [](double cutoff) { return MyFunctor(cutoff); }
   );
   ```
3. **Verify and Benchmark**:
   ```bash
   ./build/AP_Functor_Bench -v -f ATM,MyFunctor -k SoATriple --n3 both
   ```

---

## AutoPas Version Configuration

AutoPas is fetched via CMake's `FetchContent`. To benchmark against a specific branch, tag, or commit of AutoPas, update `GIT_TAG` in `cmake/modules/autopas.cmake`:

```cmake
FetchContent_Declare(
    autopasfetch
    GIT_REPOSITORY ${autopasRepoPath}
    GIT_TAG feature/3xa/atm-soa # Branch, tag, or commit hash
)
```
