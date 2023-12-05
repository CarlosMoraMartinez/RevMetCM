process trimReads{
  label 'ilm02_trimIllumina'
  conda params.trimReads.conda
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ilm02_trimIllumina", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(fastq))
  val quality
  val min_length
  
  output:
  tuple(val(illumina_id), path('*.trim.fastq.gz'))

  shell:
  '''
  trim_galore -q !{quality} --length !{min_length} --gzip --paired !{fastq[0]} !{fastq[1]}
  #rename "s/_val_[12].fq.gz/.trim.fastq.gz/" *.fq.gz #doesn't  work...
  r1=$(ls *_val_1.fq.gz)
  r2=$(ls *_val_2.fq.gz)
  mv $r1 ${r1%_val_1.fq.gz}'.trim.fastq.gz'
  mv $r2 ${r2%_val_2.fq.gz}'.trim.fastq.gz'
  '''
}