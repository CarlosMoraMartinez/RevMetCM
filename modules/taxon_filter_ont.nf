process taxonFilterOnt {
  label 'ont06_taxonFiltOnt'
  conda params.taxonFilterOnt.conda
  cpus params.resources.taxonFilterOnt.cpus
  memory params.resources.taxonFilterOnt.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont06_taxonFiltOnt", mode: 'symlink'

  input:
  path fastq
  path db
  
  output:
  tuple(path("*taxons.txt.gz"), path("*.fastq.gz"))

  shell:
  '''
  outfile=$(basename -s .fastq.gz !{fastq})
  
  kraken2 --db !{db} --threads !{params.resources.taxonFilterOnt.cpus} --unclassified-out $outfile'.tx.fastq' --gzip-compressed !{fastq} | gzip > $outfile'_taxons.txt.gz' 
  gzip $outfile'.tx.fastq' #TO DO: include pigz in the conda environment

  '''
}