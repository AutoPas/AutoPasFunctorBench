#!/bin/bash
#SBATCH -J PatternBenchmark
#SBATCH -D  ./
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=04:00:00
module load slurm_setup

export OMP_NUM_THREADS=1

./AutoPasFunctorBench

