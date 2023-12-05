process organelleGenomeIndex {
  label 'ont07_organelleFilterOnt'
  conda params.organelleFilter.conda
  cpus params.resources.organelleFilterIndex.cpus
  memory params.resources.organelleFilterIndex.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont07_organelleFilterOnt", mode: 'symlink'

  input:
  path fasta
  
  output:
  path("*.mmi")

  shell:
  '''
  outfile=$(basename -s .fna.gz !{fasta})'.mmi'
  
  minimap2 -x sr -d $outfile !{fasta}
  '''
}