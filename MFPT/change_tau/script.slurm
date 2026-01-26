#!/bin/bash 
#SBATCH --job-name=tau_MFPT_1000
#SBATCH --time=6:0:0 
#SBATCH --account=bcamley1
#SBATCH --mail-type=end 
#SBATCH --mail-user=wwang138@jh.edu 
#SBATCH --ntasks=500
#SBATCH --cpus-per-task=1 


### exit on error 
set -e 

#module load python

mpiexec python MFPT.py

wait 
date 
