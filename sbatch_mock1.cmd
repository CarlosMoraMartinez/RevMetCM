#!/bin/bash
#SBATCH -J test_job
#SBATCH -q normal
#SBATCH -D .
#SBATCH -p genD
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 128
#SBATCH --mem-per-cpu=2GB
#SBATCH -t 20:00:00
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=carlos.mora@cnag.crg.eu



nextflow revmet2.nf -c conf_crg_mock1.config -profile conda -resume
