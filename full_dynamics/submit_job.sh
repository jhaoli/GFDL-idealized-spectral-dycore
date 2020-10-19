#!/bin/bash
#SBATCH --job-name gfdl_hs
#SBATCH --partition i01203share
#SBATCH --ntasks 32
#SBATCH --nodes 2
#SBATCH --exclusive
#SBATCH --ntasks-per-node=16
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH -t 7-0:00

./run_full_dynamics_on_taiyuan.csh
