include { getIlluminaSampleList } from '../modules/get_illumina_sample_list'
include { trimReads } from '../modules/trim_reads'
include { getFastQCIllumina } from '../modules/get_fastqc_illumina'
include { taxonFilterIllumina } from '../modules/taxon_filter_illumina'


workflow ILLUMINAFASTQ {
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