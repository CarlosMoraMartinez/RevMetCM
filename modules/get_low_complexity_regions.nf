process getLowComplexityRegions {
  label 'ont05_dust'
  conda params.getLowComplexityRegions.conda
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont05_dust", mode: 'symlink'

  input:
  path fasta_file
  val window
  val level
  
  output:
  path('*.bed')

  shell:
  '''
  outfile=$(basename -s .fasta.gz !{fasta_file})
  dustmasker -in <(zcat !{fasta_file}) -outfmt acclist -level !{level} -window !{window} > $outfile'.intlist'
  paste <(awk '{print $1}' $outfile'.intlist') <(cut -f 2,3 $outfile'.intlist') | sed "s/>//" > $outfile'.bed'
  rm $outfile'.intlist'
  '''
}