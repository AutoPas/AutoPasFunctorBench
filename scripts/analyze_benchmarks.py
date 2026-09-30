#!/usr/bin/env python3
"""
AutoPas Functor Benchmark Analyzer & Plotter

Ingests Google Benchmark JSON output from AutoPasFunctorBench, parses functor/kernel/N3
metadata, normalizes performance against a single designated Baseline functor, and generates
summary tables and publication-quality plots.

Usage Examples:
    # Single JSON file (contains baseline and candidate(s)):
    python3 scripts/analyze_benchmarks.py results.json --baseline ATM --output-dir plots/

    # Multiple JSON files across optimization iterations:
    python3 scripts/analyze_benchmarks.py --baseline ATM=results_baseline.json \
        --candidate "Unroll=results_iter1.json" --candidate "AVX512=results_iter2.json" \
        --output-dir plots/
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Set styling
sns.set_theme(style="whitegrid", palette="tab10")
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
})


KNOWN_KERNELS = ["AoS", "SoASingle", "SoAPair", "SoATriple"]


def parse_benchmark_name(full_name: str) -> Dict:
    """
    Parses benchmark string like:
    BM_ATM_SoATriple_N3ON/16/3/3 or BM_ATM2_AoS/32/3/3
    """
    parts = full_name.split("/")
    bench_prefix = parts[0]

    # Clean aggregate suffix from parameters if present (e.g. "3_mean", "3_median")
    clean_parts = []
    for p in parts[1:]:
        for agg in ["_mean", "_median", "_stddev", "_cv"]:
            if p.endswith(agg):
                p = p[:-len(agg)]
                break
        clean_parts.append(p)

    num_particles = int(clean_parts[0]) if len(clean_parts) > 0 and clean_parts[0].isdigit() else None
    cell_size = float(clean_parts[1]) if len(clean_parts) > 1 else None
    cutoff = float(clean_parts[2]) if len(clean_parts) > 2 else None

    # Strip leading "BM_"
    clean_prefix = bench_prefix
    if clean_prefix.startswith("BM_"):
        clean_prefix = clean_prefix[3:]

    # Detect Newton-3 suffix
    n3_status = "ON"  # Default if not specified
    if clean_prefix.endswith("_N3ON"):
        n3_status = "ON"
        clean_prefix = clean_prefix[:-5]
    elif clean_prefix.endswith("_N3OFF") or clean_prefix.endswith("_NoN3"):
        n3_status = "OFF"
        clean_prefix = clean_prefix[:-6] if clean_prefix.endswith("_N3OFF") else clean_prefix[:-5]
    elif clean_prefix.endswith("_N3"):
        n3_status = "ON"
        clean_prefix = clean_prefix[:-3]

    # Find which known kernel matches the end of clean_prefix
    matched_kernel = None
    functor_name = clean_prefix

    for k in sorted(KNOWN_KERNELS, key=len, reverse=True):
        if clean_prefix.endswith("_" + k):
            matched_kernel = k
            functor_name = clean_prefix[:-(len(k) + 1)]
            break

    if matched_kernel is None:
        # Fallback heuristic: split by underscore
        tokens = clean_prefix.split("_")
        if len(tokens) >= 2:
            functor_name = tokens[0]
            matched_kernel = tokens[1]
        else:
            matched_kernel = "Unknown"

    return {
        "Functor": functor_name,
        "Kernel": matched_kernel,
        "Newton3": n3_status,
        "Particles": num_particles,
        "CellSize": cell_size,
        "Cutoff": cutoff,
    }


def load_benchmark_json(file_path: str, label_override: Optional[str] = None) -> pd.DataFrame:
    """Loads a single Google Benchmark JSON file into a pandas DataFrame."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Benchmark file not found: {file_path}")

    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    benchmarks = data.get("benchmarks", [])
    context = data.get("context", {})
    branch = context.get("autopas_branch", "unknown")
    commit = context.get("autopas_commit", "unknown")[:8]

    # If the file contains aggregate statistics (mean, median, stddev), only keep the 'mean'
    has_aggregates = any(b.get("run_type") == "aggregate" for b in benchmarks)

    rows = []
    for b in benchmarks:
        if has_aggregates:
            if b.get("run_type") != "aggregate" or b.get("aggregate_name") != "mean":
                continue

        name = b.get("name", "")
        parsed = parse_benchmark_name(name)

        if label_override:
            parsed["Functor"] = label_override

        row = {
            **parsed,
            "RawName": name,
            "RealTime_ns": b.get("real_time", np.nan),
            "CPUTime_ns": b.get("cpu_time", np.nan),
            "TimeUnit": b.get("time_unit", "ns"),
            "TripletsPerSec": b.get("Triplets/s", np.nan),
            "TimePerTriplet_s": b.get("Time/Triplet", np.nan),
            "TimePerInteraction_s": b.get("Time/Interaction", np.nan),
            "HitRate_pct": b.get("Hit%", np.nan),
            "Branch": branch,
            "Commit": commit,
            "SourceFile": path.name,
        }

        # Convert Time per Triplet to picoseconds for intuitive reading
        if not np.isnan(row["TimePerTriplet_s"]):
            row["TimePerTriplet_ps"] = row["TimePerTriplet_s"] * 1e12
        else:
            row["TimePerTriplet_ps"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def load_dataset(args: argparse.Namespace) -> Tuple[pd.DataFrame, str]:
    """
    Loads benchmark data from CLI arguments and resolves the baseline functor name.
    """
    dfs = []
    baseline_functor = "ATM"

    if args.files:
        for f in args.files:
            dfs.append(load_benchmark_json(f))
        if args.baseline:
            baseline_functor = args.baseline

    if args.baseline_file:
        # Syntax: name=path or path
        if "=" in args.baseline_file:
            name, path = args.baseline_file.split("=", 1)
        else:
            name, path = "Baseline", args.baseline_file
        baseline_functor = name
        dfs.append(load_benchmark_json(path, label_override=name))

    if args.candidate_files:
        for c in args.candidate_files:
            if "=" in c:
                name, path = c.split("=", 1)
            else:
                name, path = Path(c).stem, c
            dfs.append(load_benchmark_json(path, label_override=name))

    if not dfs:
        raise ValueError("No input benchmark JSON files provided. Use --help for usage.")

    df = pd.concat(dfs, ignore_index=True)

    # Validate that baseline exists in the dataset
    all_functors = df["Functor"].unique().tolist()
    if baseline_functor not in all_functors:
        # Check case-insensitive match
        for f in all_functors:
            if f.lower() == baseline_functor.lower():
                baseline_functor = f
                break
        else:
            print(f"[WARNING] Baseline '{baseline_functor}' not found in functors: {all_functors}. "
                  f"Using '{all_functors[0]}' as baseline.")
            baseline_functor = all_functors[0]

    return df, baseline_functor


def compute_speedups(df: pd.DataFrame, baseline_functor: str) -> pd.DataFrame:
    """
    Computes speedup for each benchmark relative to the baseline functor:
    speedup = baseline_real_time / candidate_real_time
    """
    df = df.copy()
    keys = ["Kernel", "Newton3", "Particles", "CellSize", "Cutoff"]

    # Filter baseline entries
    baseline_df = df[df["Functor"] == baseline_functor]
    if baseline_df.empty:
        df["Speedup"] = np.nan
        df["BaselineTime_ns"] = np.nan
        return df

    # Map keys to baseline time
    baseline_map = baseline_df.groupby(keys)["RealTime_ns"].mean().to_dict()

    def get_baseline_time(row):
        key = tuple(row[k] for k in keys)
        return baseline_map.get(key, np.nan)

    df["BaselineTime_ns"] = df.apply(get_baseline_time, axis=1)
    df["Speedup"] = df["BaselineTime_ns"] / df["RealTime_ns"]

    return df


def print_summary_table(df: pd.DataFrame, baseline_functor: str):
    """Prints a clean ASCII/Markdown summary table comparing candidates against baseline."""
    candidates = [f for f in df["Functor"].unique() if f != baseline_functor]

    print("\n" + "=" * 80)
    print(f"AutoPas Functor Benchmark Comparison (Baseline: {baseline_functor})")
    print("=" * 80)

    for kernel in df["Kernel"].unique():
        for n3 in df["Newton3"].unique():
            subset = df[(df["Kernel"] == kernel) & (df["Newton3"] == n3)]
            if subset.empty:
                continue

            print(f"\n--- Kernel: {kernel} | Newton3: {n3} ---")
            header = f"{'Particles':>10} | {baseline_functor + ' (ns)':>14} | " + " | ".join(
                [f"{c + ' (ns)':>14} | {c + ' Speedup':>12}" for c in candidates]
            )
            print(header)
            print("-" * len(header))

            particles_list = sorted(subset["Particles"].dropna().unique())
            for p in particles_list:
                p_sub = subset[subset["Particles"] == p]
                base_row = p_sub[p_sub["Functor"] == baseline_functor]
                base_time = base_row["RealTime_ns"].values[0] if not base_row.empty else np.nan

                row_str = f"{int(p):>10} | {base_time:>14.1f} | "
                cand_strs = []
                for c in candidates:
                    c_row = p_sub[p_sub["Functor"] == c]
                    if not c_row.empty:
                        c_time = c_row["RealTime_ns"].values[0]
                        c_sp = c_row["Speedup"].values[0]
                        cand_strs.append(f"{c_time:>14.1f} | {c_sp:>11.2f}x")
                    else:
                        cand_strs.append(f"{'-':>14} | {'-':>12}")
                row_str += " | ".join(cand_strs)
                print(row_str)

    print("=" * 80 + "\n")


def plot_speedup_bar(df: pd.DataFrame, baseline_functor: str, output_path: Path):
    """
    Plots a grouped bar chart of Candidate Speedup relative to Baseline (1.0x line).
    Only candidate functors are plotted as bars; Baseline is represented by the 1.0x line.
    """
    candidates_df = df[df["Functor"] != baseline_functor].copy()
    if candidates_df.empty:
        print("[INFO] No candidate functors found to plot speedup against baseline.")
        return

    kernels = candidates_df["Kernel"].unique()
    num_kernels = len(kernels)

    fig, axes = plt.subplots(
        1, num_kernels, figsize=(6 * num_kernels, 5), squeeze=False, sharey=True
    )

    for idx, kernel in enumerate(kernels):
        ax = axes[0, idx]
        k_df = candidates_df[candidates_df["Kernel"] == kernel].copy()

        # Combine Functor + Newton3 for hue if both N3 modes exist
        if len(k_df["Newton3"].unique()) > 1:
            k_df["Variant"] = k_df["Functor"] + " (N3 " + k_df["Newton3"] + ")"
        else:
            k_df["Variant"] = k_df["Functor"]

        sns.barplot(
            data=k_df,
            x="Particles",
            y="Speedup",
            hue="Variant",
            ax=ax,
            edgecolor="black",
            linewidth=0.8,
        )

        ax.axhline(1.0, color="crimson", linestyle="--", linewidth=1.5, label=f"Baseline ({baseline_functor})")
        ax.set_title(f"Kernel: {kernel}", fontweight="bold")
        ax.set_xlabel("Particles per Cell")
        if idx == 0:
            ax.set_ylabel("Speedup vs. Baseline (Higher is Better)")
        ax.legend(title=None, loc="best")
        ax.grid(axis="y", linestyle=":", alpha=0.7)

    plt.suptitle(f"Candidate Speedup Normalized to Baseline ({baseline_functor} = 1.0x)", y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] Speedup bar plot: {output_path}")


def plot_throughput_scaling(df: pd.DataFrame, output_path: Path):
    """
    Plots throughput (GigaTriplets/sec) vs. Particle Count per Cell on a log-x scale.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    plot_df = df.copy()
    plot_df["GTripletsPerSec"] = plot_df["TripletsPerSec"] / 1e9
    plot_df["Variant"] = plot_df["Functor"] + " (" + plot_df["Kernel"] + ", N3 " + plot_df["Newton3"] + ")"

    sns.lineplot(
        data=plot_df,
        x="Particles",
        y="GTripletsPerSec",
        hue="Variant",
        style="Variant",
        markers=True,
        dashes=False,
        linewidth=2,
        markersize=8,
        ax=ax,
    )

    ax.set_xscale("log", base=2)
    ax.set_xlabel("Particles per Cell (Log Scale)")
    ax.set_ylabel("Throughput (GigaTriplets / s)")
    ax.set_title("Functor Throughput Scaling", fontweight="bold")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(title=None, bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] Throughput scaling plot: {output_path}")


def plot_time_per_triplet(df: pd.DataFrame, output_path: Path):
    """
    Plots Time per Triplet (in picoseconds) vs. Particle Count.
    Lower is better.
    """
    valid_df = df.dropna(subset=["TimePerTriplet_ps"]).copy()
    if valid_df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    valid_df["Variant"] = valid_df["Functor"] + " (" + valid_df["Kernel"] + ", N3 " + valid_df["Newton3"] + ")"

    sns.lineplot(
        data=valid_df,
        x="Particles",
        y="TimePerTriplet_ps",
        hue="Variant",
        style="Variant",
        markers=True,
        dashes=False,
        linewidth=2,
        markersize=8,
        ax=ax,
    )

    ax.set_xlabel("Particles per Cell")
    ax.set_ylabel("Time per Triplet (ps, Lower is Better)")
    ax.set_title("Time per Triplet vs. Particle Count", fontweight="bold")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(title=None, bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] Time per triplet plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="AutoPas Functor Benchmark Analyzer & Plotter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="One or more Google Benchmark JSON files to analyze.",
    )
    parser.add_argument(
        "--baseline",
        default="ATM",
        help="Name of the baseline functor (default: 'ATM'). Everything will be normalized against this functor.",
    )
    parser.add_argument(
        "--baseline-file",
        help="Specify a dedicated JSON file for the baseline (e.g. --baseline-file ATM=baseline.json).",
    )
    parser.add_argument(
        "--candidate",
        dest="candidate_files",
        action="append",
        help="Specify candidate JSON files (e.g. --candidate 'Unroll=iter1.json'). Can be used multiple times.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="plots",
        help="Directory to save generated plot images (default: 'plots/').",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Only print the summary table to stdout without generating plot files.",
    )

    args = parser.parse_args()

    # Load data
    try:
        df, baseline = load_dataset(args)
    except Exception as e:
        print(f"[ERROR] Failed to load benchmarks: {e}", file=sys.stderr)
        sys.exit(1)

    # Compute speedup relative to single baseline
    df = compute_speedups(df, baseline)

    # Print summary table
    print_summary_table(df, baseline)

    # Generate plots if requested
    if not args.no_plots:
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        plot_speedup_bar(df, baseline, out_dir / "speedup_vs_baseline.png")
        plot_throughput_scaling(df, out_dir / "throughput_scaling.png")
        plot_time_per_triplet(df, out_dir / "time_per_triplet.png")


if __name__ == "__main__":
    main()
