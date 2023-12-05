process getFastaIDs {
  label 'ont03_fids'
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  publishDir "$results_dir/ont03_getFastaIDs", mode: 'symlink'
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  
  input:
  path fasta_file
  
  output:
   path('*.ids')


  shell:
  '''
  #Get list of nanopore read ids from nanopore fasta
  outfile=$(basename -s .fasta.gz !{fasta_file})
  zgrep ">" !{fasta_file}  | sed 's/>//g' | cut -f 1 -d\\  > $outfile'.ids'
  '''
}