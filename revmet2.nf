#!/usr/bin/env nextflow

  ch_ont = Channel
        .fromPath(params.ont)
        //.view{"Input ONT: $it"}

  ch_illumina = Channel
        .fromFilePairs(params.illumina)
        //.view{"Input Illumina: $it"}

process getIlluminaSampleList{
    label 'nf_01_ilslst'
    cpus params.resources.standard1.cpus
    memory params.resources.standard1.mem
    errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'terminate' }
    maxRetries 10
    publishDir "$results_dir/1_getIlluminaSampleList"
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

process filterOntReads {
    label 'nf_02_flng'
    conda params.filterOntReads.conda
    cpus params.resources.standard2.cpus
    memory params.resources.standard2.mem
    errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
    maxRetries 10
    publishDir "$results_dir/2_filterOntReads-l$min_length-q$min_mean_q", mode: 'copy'
    input:
    path ont_file
    val min_length
    val min_mean_q
    
    output:
    path('*.fastq.gz')

    shell:
    '''
    #Get fasta from nanopore fastq
    outfile=$(basename -s .fastq.gz !{ont_file} | sed "s/_/-/g")
    filtlong --min_length !{min_length} --min_mean_q !{min_mean_q} !{ont_file} | gzip > $outfile'.trim.fastq.gz'

    '''

}

process makeFastaFromFastq {
    label 'nf_03_gtfa'
    cpus params.resources.standard1.cpus
    memory params.resources.standard1.mem
    errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
    maxRetries 10
    publishDir "$results_dir/3_makeFastaFromFastq", mode: 'symlink'
    input:
    path ont_file
    
    output:
    path('*.fasta.gz')

    shell:
    '''
    #Get fasta from nanopore fastq
    outfile=$(basename -s .fastq.gz !{ont_file} | sed "s/_/-/g")
    seqtk seq -a !{ont_file} | gzip > $outfile'.fasta.gz'
    '''
}

process getFastaIDs {
    label 'nf_04_fids'
    cpus params.resources.standard1.cpus
    memory params.resources.standard1.mem
    publishDir "$results_dir/4_getFastaIDs", mode: 'symlink'
    errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
    maxRetries 10
    input:
    path fasta_file
    
    output:
     path('*.ids')


    shell:
    '''
    #Get list of nanopore read ids from nanopore fasta
    outfile=$(basename -s .fasta.gz !{fasta_file})
    zgrep ">" !{fasta_file}  | sed 's/>//g' | cut -f 1 -d\\  > $outfile'.ids'
    '''

}

process indexReferenceBwa {
    label 'nf_05_idx'
    cpus params.resources.index.cpus
    memory params.resources.index.mem
    errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
    maxRetries 10
    publishDir "$results_dir/5_ont_index_bwa", mode: 'symlink'

    input:
    path fasta_file
    
    output:
    tuple path("${fasta_file}"), path("${fasta_file}*")


    shell:
    '''
    #Index nanopore fasta reference prior to alignment
    #Note: this will happen even if minimap2 is used instead of bwa. Fix in the future.
    bwa index !{fasta_file}
    '''
}

process alignIllumina {
  label 'nf_06_algn'
  cpus params.resources.alignment.cpus
  memory params.resources.alignment.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'terminate' }
  maxRetries 10
  publishDir "$results_dir/6_alignIllumina", mode: 'symlink'
  input:
    tuple(val(ont_file), val(ont_index), val(illumina_id), val(illumina_reads))
    
  output:
    path('*_raw.bam')
    
  shell:
  if( params.alignIllumina.program == 'minimap2' )
      '''
        outfile=$(basename -s .fasta !{ont_file})_!{illumina_id}
        rawbam=$outfile'_raw.bam'
  
        minimap2 -ax sr !{ont_file} !{illumina_reads[0]}  !{illumina_reads[1]} | samtools view -bS -o $rawbam
      '''
  else if( params.alignIllumina.program == 'bwa' )
      '''
        outfile=$(basename -s .fasta.gz !{ont_file})_!{illumina_id}
        rawbam=$outfile'_raw.bam'

        bwa mem -t !{params.resources.alignment.cpus} !{params.alignIllumina.bwaparams} !{ont_file} !{illumina_reads[0]}  !{illumina_reads[1]} | samtools view -bS -o $rawbam
      '''
}

process filterIlluminaAlignment{
  label 'nf_07_smflt'
  cpus params.resources.samtoolsfilter.cpus
  memory params.resources.samtoolsfilter.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/7_filterIlluminaAlignment-F$exclude_flag_F-q$mapq", mode: 'symlink'
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

process getCoverage{
  label 'nf_08_dpth'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'terminate' }
  maxRetries 10
  publishDir "$results_dir/8_getCoverage", mode: 'symlink'
  input:
    path bam

  output:
    path('*.depth.gz')

  shell:
  '''
      outfile=$(basename -s .bam !{bam})
      samtoolsdepth=$outfile'.depth.gz'

      #Get coverage at each position on each nanopore read
      samtools coverage --ff 0 !{bam} | gzip > $samtoolsdepth
  '''
}

process mergeAndBin2species {
  label 'nf_09_bin2sp'
  conda params.mergeAndBin2species.conda
  cpus params.resources.mergeandbin2species.cpus
  memory params.resources.mergeandbin2species.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  storeDir "$results_dir/9_MergeAndBin"
  input:
    tuple(val(ont_id), path(pcs))

  output:
    path "${ont_id}_all.csv"

  shell:
  '''
      python !{params.scriptsdir}/merge_coverages.py  -n !{params.resources.standard2.cpus} -p !{ont_id}
  '''
}

workflow {

  if(params.filterOntReads.do_filter){
    filterOntReads(ch_ont, params.filterOntReads.min_length, params.filterOntReads.min_mean_q) |
    makeFastaFromFastq
  }else{
    makeFastaFromFastq(ch_ont)
  }
  ch_ont_fastq = makeFastaFromFastq.out
  //ch_ont_fastq.view{ "Fasta created from fastq: $it" }

  indexReferenceBwa(ch_ont_fastq)
  ch_ont_index = indexReferenceBwa.out
  //ch_ont_index.view{ "BWA index created from fasta: $it" }

  getFastaIDs(ch_ont_fastq)
  ch_ont_ids = getFastaIDs.out
  //ch_ont_ids.view{ "Fasta ID list created from fasta: $it" }

  ch_illumina_samples = ch_illumina
    .map{x ->
      def sname = x.get(0)
      return sname}
    .collect()
    //.view()
  getIlluminaSampleList(ch_illumina_samples)
  ch_illumina_samplelist=getIlluminaSampleList.out
  //.view{ "Illumina sample list: $it" }

  ch_ont_index=ch_ont_index.combine(ch_illumina)
    //.view{ "ch_ont_index combined: $it" }
  
  alignIllumina(ch_ont_index) 
  ch_aligned = alignIllumina.out
    //.view{ "Alignment result: $it" }
  filterIlluminaAlignment(
    ch_aligned,
    params.filterIlluminaAlignment.mapq,
    params.filterIlluminaAlignment.include_flag_f,
    params.filterIlluminaAlignment.exclude_flag_F
  ) | 
  //view{ "filterIlluminaAlignmentignment result: $it" } |
  getCoverage 
  //| view{ "getCoverage result: $it" } 

  ch_pcs = getCoverage.out
    .map{ file ->
      def key = file.name.toString().tokenize('_').get(0)
      return tuple(key, file)
    }
    .groupTuple()
  .view{ "getCoverage grouped by nanopore sample: $it" }

  //Merge all pcs from one nanopore sample into one file
  mergeAndBin2species(ch_pcs)
    .view{ "mergeAndBin2species results: $it" }

}