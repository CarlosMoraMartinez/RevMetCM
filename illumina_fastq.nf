process getIlluminaSampleList{
  label 'ilm01_ilslst'
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ilm01_getIlluminaSampleList"
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

process taxonFilterIllumina {
  label 'ilm04_taxonFiltIllumina'
  conda params.taxonFilterIllumina.conda
  cpus params.resources.taxonFilterIllumina.cpus
  memory params.resources.taxonFilterIllumina.mem
  errorStrategy { task.exitStatus in 1..2 ? 'terminate' : 'terminate' }
  maxRetries 0
  publishDir "$results_dir/ilm04_taxonFiltIllumina", mode: 'symlink'

  input:
  tuple(val(illumina_id), path(fastq))
  path db
  
  output:
  tuple(path("*taxons.txt.gz"), val(illumina_id), path('*.tx.fastq.gz'))

  shell:
  ''' 
  outfile=!{illumina_id}'_taxons.txt.gz' 
  kraken2 --db !{db} --paired --threads !{params.resources.taxonFilterOnt.cpus} --unclassified-out unclassified#.fastq --gzip-compressed !{fastq[0]} !{fastq[1]} | gzip > $outfile
  gzip unclassified_1.fastq unclassified_2.fastq
  mv unclassified_1.fastq.gz $(basename -s .fastq.gz !{fastq[0]}).tx.fastq.gz
  mv unclassified_2.fastq.gz $(basename -s .fastq.gz !{fastq[1]}).tx.fastq.gz
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
    //.view{ "Illumina trimmed reads: $it" }

  ch_taxons_illumina = Channel.from([])
  if(params.taxonFilterIllumina.do_filter){
    taxonFilterIllumina(ch_illumina_processed, params.taxonFilterIllumina.db)
    ch_illumina_processed = taxonFilterIllumina.out.map{it -> [it[1], it[2]]}
        //.view{ "Illumina fastq filtered by taxon: $it" }
    ch_taxons_illumina = taxonFilterIllumina.out.map{it -> it[0]}
        //.view{ "Illumina taxon composition: $it" }
  }

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
    ch_taxons_illumina

}