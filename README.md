# Reverse metagenomics pipeline to quantify species in pollen samples. 

Briefly, Illumina reads from individual species are mapped against long reads from environmental samples, and each of these are assigned to the species that covers it the most. 

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

