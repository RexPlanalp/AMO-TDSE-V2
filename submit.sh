#!/usr/bin/bash
#SBATCH --job-name testing
#SBATCH --output run.log
#SBATCH --nodes 1
#SBATCH --ntasks 32
#SBATCH --mem=8G

#SBATCH --exclude=node83
#SBATCH --exclude=node81

source ~/.bashrc

REPO_DIR="/users/becker/dopl4670/Research/AMO-TDSE/build/bin"
hostname
pwd

chmod +x $REPO_DIR/simulation.exe

mpirun -np $SLURM_NTASKS $REPO_DIR/simulation.exe $@ >> results.log

