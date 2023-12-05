process alignIllumina {
  label 'rvm01_algn'
  cpus params.resources.alignment.cpus
  memory params.resources.alignment.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/rvm01_alignIllumina", mode: 'symlink'

  input:
    tuple(val(ont_file), val(ont_index), val(illumina_id), val(illumina_reads))

  output:
    path('*_raw.bam')

  shell:
  if( params.alignIllumina.program == 'minimap2' )
    '''
    #outfile=$(basename -s .fasta.gz !{ont_file})_!{illumina_id}
    outfile=$(cut -f 1 -d. <(basename !{ont_file}))_!{illumina_id}
    rawbam=$outfile'_raw.bam'

    minimap2 -t !{params.resources.alignment.cpus} -ax sr !{ont_index} !{illumina_reads[0]}  !{illumina_reads[1]} | samtools view -bS -o $rawbam
    '''
  else if( params.alignIllumina.program == 'bwa' )
    '''
    outfile=$(cut -f 1 -d. <(basename !{ont_file}))_!{illumina_id}
    rawbam=$outfile'_raw.bam'

    bwa mem -t !{params.resources.alignment.cpus} !{params.alignIllumina.bwaparams} !{ont_file} !{illumina_reads[0]}  !{illumina_reads[1]} | samtools view -bS -o $rawbam
    '''
}