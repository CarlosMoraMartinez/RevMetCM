#!/bin/bash
#SBATCH -J mrgpy
#SBATCH -q short
#SBATCH -D .
#SBATCH -p genD
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem-per-cpu=2GB
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=carlos.mora@cnag.crg.eu

echo $1
python merge_coverages.py -m 15 -n 100 -p $1 -o binresults
