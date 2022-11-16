#!/bin/bash
#SBATCH -J mock1
#SBATCH -q normal
#SBATCH -D .
#SBATCH -p genD
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 2
#SBATCH --mem-per-cpu=2GB
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=carlos.mora@cnag.crg.eu



nextflow revmet2.nf -c FAR77735_conf_crg_mock1.config -profile conda -resume
