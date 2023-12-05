process organelleFilterOnt {
  label 'ont08_organelleFilterOnt'
  conda params.organelleFilter.conda
  cpus params.resources.organelleFilter.cpus
  memory params.resources.organelleFilter.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont08_organelleFilterOnt", mode: 'symlink'

  input:
  path index
  path fastq
  
  output:
  path("*.nuc.fastq.gz")

  shell:
  '''
  outfile=$(basename -s .fastq.gz !{fastq})'.nuc.fastq.gz'
  
  minimap2 -t !{params.resources.organelleFilter.cpus} -ax sr !{index} !{fastq} | samtools view -f 4 | cut -f 1 > ids.txt
  seqtk subseq !{fastq} ids.txt | gzip > $outfile

  '''
}