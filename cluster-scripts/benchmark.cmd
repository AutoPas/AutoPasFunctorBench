#!/bin/bash
#SBATCH -J functor_benchmark
#SBATCH -o ./output/functor_bench_%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mail-type=end
#SBATCH --mail-user=luis.gall@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:05:00
module load slurm_setup
module load gcc
module load cmake

./build/AutoPasFunctorBench