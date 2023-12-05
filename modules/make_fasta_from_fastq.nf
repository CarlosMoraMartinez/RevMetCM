process makeFastaFromFastq {
  label 'ont02_gtfa'
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont02_makeFastaFromFastq", mode: 'symlink'
  
  input:
  path ont_file
  
  output:
  path('*.fasta.gz')

  shell:
  '''
  #Get fasta from nanopore fastq
  outfile=$(basename -s .fastq.gz !{ont_file} | sed "s/_/-/g")
  seqtk seq -a !{ont_file} | gzip > $outfile'.fasta.gz'
  '''
}