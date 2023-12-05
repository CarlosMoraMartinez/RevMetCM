process getFastQCIllumina{
  label 'ilm03_fastqc'
  conda params.getFastQCIllumina.conda
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ilm03_FastQC", mode: 'copy'
  input:
  path fastq
  
  output:
  path '*.html'

  shell:
  '''
  fastqc -q !{fastq}  
  '''
}