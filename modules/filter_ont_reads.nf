process filterOntReads {
  label 'ont01_flng'
  conda params.filterOntReads.conda
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont01_filterOntReads-l$min_length-q$min_mean_q", mode: 'copy'
  
  input:
  path ont_file
  val min_length
  val min_mean_q
  
  output:
  path('*trim.fastq.gz')

  shell:
  '''
  #Get fasta from nanopore fastq
  outfile=$(basename -s .fastq.gz !{ont_file} | sed "s/_/-/g")
  filtlong --min_length !{min_length} --min_mean_q !{min_mean_q} !{ont_file} | gzip > $outfile'.trim.fastq.gz'
  '''
}