#!/usr/bin/env nextflow
include { ont2fasta } from './process_ont.nf'
include { illuminafastq } from './illumina_fastq.nf'

ch_ont = Channel
  .fromPath(params.ont)
  //.view{"Input ONT: $it"}

ch_illumina = Channel
  .fromFilePairs(params.illumina)
  //.view{"Input Illumina: $it"}


process alignIllumina {
  label 'nf_09_algn'
  cpus params.resources.alignment.cpus
  memory params.resources.alignment.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/09_alignIllumina", mode: 'symlink'

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

process filterIlluminaAlignment{
  label 'nf_10_smflt'
  cpus params.resources.samtoolsfilter.cpus
  memory params.resources.samtoolsfilter.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/10_filterIlluminaAlignment-F$exclude_flag_F-q$mapq", mode: 'symlink'

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

process filterDustRegions{
  label 'nf_11_dustfilt'
  conda params.filterDustRegions.conda
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/11_filterDust", mode: 'symlink'
  
  input:
    tuple(val(ont), val(bam), val(bed))

  output:
    path('*.dust.bam')

  shell:
  '''
  outfile=$(basename -s .bam !{bam})
  outbam=$outfile'.dust.bam'

  bedtools intersect -v -abam !{bam} -b !{bed} > $outbam
  '''
}

process getCoverage{
  label 'nf_12_dpth'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/12_getCoverage", mode: 'symlink'
  
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
  label 'nf_13_bin2sp'
  conda params.mergeAndBin2species.conda
  cpus params.resources.mergeandbin2species.cpus
  memory params.resources.mergeandbin2species.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 3
  storeDir "$results_dir/13_MergeAndBin"

  input:
    tuple(val(ont_id), path(pcs))
    val min_coverage

  output:
    tuple(path("${ont_id}*_all.csv"), path("${ont_id}_tmp.txt"))

  shell:
  '''
  for pc in !{pcs}; do echo $pc >> !{ont_id}_tmp.txt; done
  python !{params.scriptsdir}/merge_coverages.py  -m !{min_coverage} -n !{params.resources.mergeandbin2species.cpus} -p !{ont_id}
  '''
}

workflow {


  ont2fasta(ch_ont)
  ch_ont_index = ont2fasta.out.ch_ont_index

  illuminafastq(ch_illumina)
  ch_illumina_samplelist = illuminafastq.out.ch_illumina_samplelist
  ch_illumina_processed = illuminafastq.out.ch_illumina_processed

  ch_ont_index=ch_ont_index.combine(ch_illumina_processed)
  //.view{ "ch_ont_index combined: $it" }
  
  alignIllumina(ch_ont_index) 
  ch_aligned = alignIllumina.out
  //.view{ "Alignment result: $it" }
  filterIlluminaAlignment(
    ch_aligned,
    params.filterIlluminaAlignment.mapq,
    params.filterIlluminaAlignment.include_flag_f,
    params.filterIlluminaAlignment.exclude_flag_F
  )
  //.view{ "filterIlluminaAlignmentignment result: $it" } |
  ch_filtered = filterIlluminaAlignment.out

  //Filter aligmnents in low complexity regions detected by dustmasker
  if(params.filterLowComplexityRegions){
    ch_dustbed = ont2fasta.out.ch_ont_dust
    .map{ file ->
      def key = file.name.toString().tokenize('.').get(0)
      return tuple(key, file)
    }
     //.view{ "filterLowComplexityRegions - map ont name to bed file: $it" } 
    ch_filtered = ch_filtered.map{ file ->
      def key = file.name.toString().tokenize('_').get(0)
      return tuple(key, file)
    }.combine(ch_dustbed, by:0)
     //.view{ "filterLowComplexityRegions - map bed with bam: $it" } 
    filterDustRegions(ch_filtered)
    ch_filtered = filterDustRegions.out
     //.view{ "filterDustRegions output bam: $it" }
  }

  getCoverage(ch_filtered)
  //| view{ "getCoverage result: $it" } 

  ch_pcs = getCoverage.out
  .map{ file ->
    def key = file.name.toString().tokenize('_').get(0)
    return tuple(key, file)
  }
  .groupTuple()
  //.view{ "getCoverage grouped by nanopore sample: $it" }

  //Merge all pcs from one nanopore sample into one file
  mergeAndBin2species(ch_pcs, params.mergeAndBin2species.min_coverage)
  //.view{ "mergeAndBin2species results: $it" }

}