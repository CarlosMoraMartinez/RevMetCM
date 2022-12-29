conda activate /software/assembly/conda/revmet/
for i in $(ls | cut -f 1 -d_ | sort | uniq); do sbatch sbatch_merge_all.cmd $i; done

