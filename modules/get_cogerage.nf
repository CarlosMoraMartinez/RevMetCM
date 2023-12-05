process getCoverage{
  label 'rvm04_dpth'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/rvm04_getCoverage", mode: 'symlink'
  
  input:
    path bam

  output:
    path('*.depth.gz')

  shell:
  '''
  outfile=$(basename -s .bam !{bam})
  samtoolsdepth=$outfile'.depth.gz'

  #Get coverage at each position on each nanopore read
  samtools coverage --ff 0 !{bam} | gzip > $samtoolsdepth
  '''
}