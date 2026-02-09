#!/bin/bash
#SBATCH --job-name=2d_AOUPs
#SBATCH --time=6:0:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=begin,end,fail,time_limit
#SBATCH --mail-user=wwang138@jh.edu
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1

### exit on error
set -e

python --version
mpiexec python scan_va.py

wait
date

