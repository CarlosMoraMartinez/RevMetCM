#!/usr/bin/env nextflow
include { ONT2FASTA } from './workflows/process_ont.nf'
include { ILLUMINAFASTQ } from './workflows/illumina_fastq.nf'
include { REVMETCORE } from './workflows/revmet_core.nf'

workflow {

  ch_ont = Channel
    .fromPath(params.ont)
    //.view{"Input ONT: $it"}

  ch_illumina = Channel
    .fromFilePairs(params.illumina)
    //.view{"Input Illumina: $it"}

  ONT2FASTA(ch_ont)
  ch_ont_index = ONT2FASTA.out.ch_ont_index
  ch_ont_dust = ONT2FASTA.out.ch_ont_dust

  ILLUMINAFASTQ(ch_illumina)
  ch_illumina_samplelist = ILLUMINAFASTQ.out.ch_illumina_samplelist
  ch_illumina_processed = ILLUMINAFASTQ.out.ch_illumina_processed

  REVMETCORE(ch_ont, ch_ont_dust, ch_illumina_processed)
  ch_merged_results = REVMETCORE.out.ch_merged
    .view{"Final pipeline result: $it"}
  
}