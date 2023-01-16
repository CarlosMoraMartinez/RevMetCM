# Reverse metagenomics pipeline to quantify species in pollen samples. 

Briefly, Illumina reads from individual species are mapped against long reads from environmental samples, and each of these are assigned to the species that covers it the most. 

Uses Kraken2 and alignment to a custom fasta database to filter out contaminant/confounding reads.

Based on the original RevMet pipeline from https://github.com/nedpeel/RevMet 

## Usage

To run locally use: 

```
nextflow run revmet2.nf -c config/config_test.config  -profile conda -resume -with-dag flowchart4.png
```

To run in a slurm cluster run the following command with the NextFlow environment active: 

```
sbatch sbatch/sbatch_samples_nuc.cmd
```

In case you want to do the read assignment only, go to the directory where the coverage files are, copy the launch_mergepy.sbatch and  and execute:

```
cp ../RevMetCM/sbatch/launch_mergepy.sbatch .
cp ../RevMetCM/scripts/merge_coverages.py .

for i in $(ls | cut -f 1 -d_ | sort | uniq); do sbatch launch_mergepy.sbatch $i; done
```

The downstream analysis was carried out locally with the *scripts/analyze_from_merged.R* R script.


## Input

Input directories are indicated in the .config file.

- Illumina reads from independent species
- ONT reads from environmental mixtures

## Output
- *.csv* files with each ONT read assigned to one of the species with short reads. 


## Software versions used

- NextFlow	22.04.0
- samtools	1.16.1
- bedtools	2.30.0
- filtlong	0.2.1
- dustmasker	1.0.0
- Kraken	2.1.2
- seqtk	1.3-r106
- Trim Galore	0.6.7
- FastQC	v0.11.9
- bwa mem	0.7.17-r1188
- Minimap2 2.17-r941
- MASH	2.3
- python	3.7.12
- pandas	1.3.5

