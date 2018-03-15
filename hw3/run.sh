#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:10:00
#SBATCH -p standard
#SBATCH --output=output.txt

# Run program
module load intel
cd $SLURM_SUBMIT_DIR
./mse627