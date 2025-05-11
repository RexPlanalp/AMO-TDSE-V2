#!/usr/bin/bash
#SBATCH --job-name testing
#SBATCH --output run.log
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --mem=8G

source ~/.bashrc

REPO_DIR="/users/becker/dopl4670/Research/AMO-TDSE/build/bin"
hostname
pwd

chmod +x $REPO_DIR/simulation.exe

mpirun -np $SLURM_NTASKS $REPO_DIR/simulation.exe $@ >> results.log

