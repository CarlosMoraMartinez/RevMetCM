process filterDustRegions{
  label 'rvm03_dustfilt'
  conda params.filterDustRegions.conda
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/rvm03_filterDust", mode: 'symlink'
  
  input:
    tuple(val(ont), val(bam), val(bed))

  output:
    path('*.dust.bam')

  shell:
  '''
  outfile=$(basename -s .bam !{bam})
  outbam=$outfile'.dust.bam'

  bedtools intersect -v -abam !{bam} -b !{bed} > $outbam
  '''
}
