#!/bin/bash 
#SBATCH --job-name=growth_models
#SBATCH --time=1:0:0 
#SBATCH --account=bcamley1
#SBATCH --mail-type=end 
#SBATCH --mail-user=wwang138@jh.edu 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 

#SBATCH -p shared

### exit on error 
set -e 

#module load python

python all_growth.py

wait 
date 