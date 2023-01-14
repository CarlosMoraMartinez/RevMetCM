#!/bin/bash
#SBATCH -J msh
#SBATCH -q short
#SBATCH -D .
#SBATCH -p genD
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem-per-cpu=2GB
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=carlos.mora@cnag.crg.eu

python sketch_all.py

