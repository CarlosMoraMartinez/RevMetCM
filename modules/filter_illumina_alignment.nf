process filterIlluminaAlignment{
  label 'rvm02_smflt'
  cpus params.resources.samtoolsfilter.cpus
  memory params.resources.samtoolsfilter.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/rvm02_filterIlluminaAlignment-F$exclude_flag_F-q$mapq", mode: 'symlink'

  input:
    path bam
    val mapq
    val include_flag_f
    val exclude_flag_F

  output:
    path('*.filt.bam')

  shell:
  '''
  outfile=$(basename -s _raw.bam !{bam})
  filtbam=$outfile'.filt.bam'

  #Filter and sort sam alignment file then output as bam
  samtools view -@ !{params.resources.samtoolsfilter.cpus} -b -F !{exclude_flag_F} -f !{include_flag_f} -q !{mapq} !{bam} | samtools sort -o $filtbam
  samtools index $filtbam
  '''
}