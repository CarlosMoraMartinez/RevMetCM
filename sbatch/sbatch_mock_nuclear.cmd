#!/bin/bash
#SBATCH -J nuc
#SBATCH -q eternal
#SBATCH -D .
#SBATCH -p genD
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -c 2
#SBATCH --mem-per-cpu=2GB
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=carmoma9@gmail.com

nextflow run revmet2.nf -c config/run_mock_cluster.config -profile conda -resume -with-report report.html
