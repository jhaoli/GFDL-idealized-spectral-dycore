#!/bin/bash
#SBATCH --job-name gfdl_test
#SBATCH --partition i01203share
#SBATCH --ntasks 2
#SBATCH --nodes 1
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH -t 7-0:00

./run_full_dynamics_on_taiyuan_rh4.csh
