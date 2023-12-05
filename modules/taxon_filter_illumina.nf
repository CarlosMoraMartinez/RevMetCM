process taxonFilterIllumina {
  label 'ilm04_taxonFiltIllumina'
  conda params.taxonFilterIllumina.conda
  cpus params.resources.taxonFilterIllumina.cpus
  memory params.resources.taxonFilterIllumina.mem
  errorStrategy { task.exitStatus in 1..2 ? 'terminate' : 'terminate' }
  maxRetries 0
  publishDir "$results_dir/ilm04_taxonFiltIllumina", mode: 'symlink'

  input:
  tuple(val(illumina_id), path(fastq))
  path db
  
  output:
  tuple(path("*taxons.txt.gz"), val(illumina_id), path('*.tx.fastq.gz'))

  shell:
  ''' 
  outfile=!{illumina_id}'_taxons.txt.gz' 
  kraken2 --db !{db} --paired --threads !{params.resources.taxonFilterOnt.cpus} --unclassified-out unclassified#.fastq --gzip-compressed !{fastq[0]} !{fastq[1]} | gzip > $outfile
  gzip unclassified_1.fastq unclassified_2.fastq
  mv unclassified_1.fastq.gz $(basename -s .fastq.gz !{fastq[0]}).tx.fastq.gz
  mv unclassified_2.fastq.gz $(basename -s .fastq.gz !{fastq[1]}).tx.fastq.gz
  '''
}