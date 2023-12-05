process indexReference {
  label 'ont04_idx'
  conda params.indexReference.conda
  cpus params.resources.index.cpus
  memory params.resources.index.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont04_index", mode: 'symlink'

  input:
  path fasta_file
  
  output:
  tuple path("${fasta_file}"), path("${fasta_file}*")


  shell:
  if(params.alignIllumina.program == "minimap2")
    '''
    minimap2 -x sr -d !{fasta_file}'.mmi' !{fasta_file}
    '''
  else if(params.alignIllumina.program == "bwa")
    '''
    #Index nanopore fasta reference prior to alignment
    #Note: this will happen even if minimap2 is used instead of bwa. Fix in the future.
    bwa index !{fasta_file}
  '''
}