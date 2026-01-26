#!/bin/bash 
#SBATCH --job-name=NAOUPs
#SBATCH --time=3:0:0 
#SBATCH --account=bcamley1
#SBATCH --mail-type=end 
#SBATCH --mail-user=wwang138@jh.edu 
#SBATCH --ntasks=192
#SBATCH --cpus-per-task=1 


### exit on error 
set -e 

#module load python

mpiexec python NAOUPs.py
#mpiexec python NAOUPs_smaller_dt.py

wait 
date
