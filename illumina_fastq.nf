process getIlluminaSampleList{
  label 'nf_01_ilslst'
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/01_getIlluminaSampleList"
  input:
  val ill_names
  
  output:
  path 'skim_ref.ids'

  shell:
  '''
  touch skim_ref.ids
  echo !{ill_names} | sed 's/, /\\n/g' | tr -d [] | cut -f 1 -d_ > skim_ref.ids
  '''
}

process trimReads{
  label 'nf_02_trimIllumina'
  conda params.trimReads.conda
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/02_trimIllumina", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(fastq))
  val quality
  val min_length
  
  output:
  tuple(val(illumina_id), path('*.trim.fastq.gz'))

  shell:
  '''
  trim_galore -q !{quality} --length !{min_length} --gzip --paired !{fastq[0]} !{fastq[1]}
  rename "s/_val_[12].fq.gz/.trim.fastq.gz/" *.fq.gz
  '''
}

process getFastQCIllumina{
  label 'nf_03_fastqc'
  conda params.getFastQCIllumina.conda
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  storeDir "$results_dir/03_FastQC"
  input:
  path fastq
  
  output:
  path '*.html'

  shell:
  '''
  fastqc -q !{fastq}  
  '''
}

workflow illuminafastq {
  take: ch_illumina
  main:
  ch_illumina_samples = ch_illumina
  .map{x ->
    def sname = x.get(0)
    return sname}
  .collect()
  //.view()
  getIlluminaSampleList(ch_illumina_samples)
  ch_illumina_samplelist=getIlluminaSampleList.out
  //.view{ "Illumina sample list: $it" }

  //Add raw fastq to a channel to perform FastQC
  ch_fastqc = Channel.from([])
  ch_fastq = Channel.from([])
  if(params.getFastQCIllumina.do_fastqc){
    ch_fastq = ch_illumina.map{it -> it[1]}.flatten()
  }

  //Trim reads
  ch_illumina_processed = ch_illumina
  if(params.trimReads.do_trim){
    trimReads(ch_illumina, params.trimReads.quality, params.trimReads.min_length)
    ch_illumina_processed = trimReads.out
    .view{ "Illumina trimmed reads: $it" }

  //Add trimmed fastq to a channel to perform FastQC
    if(params.getFastQCIllumina.do_fastqc){
      ch_fastq = ch_illumina_processed.map{it -> it[1]}
      .flatten().concat(ch_fastq)
    }
  }

  //Perform FastQC
  if(params.getFastQCIllumina.do_fastqc){
    getFastQCIllumina(ch_fastq)
    //.view{ "Illumina FastQC reports: $it" }
    ch_fastqc = getFastQCIllumina.out
  }

  emit:
    ch_illumina_samplelist
    ch_fastqc
    ch_illumina_processed

}