#!/usr/bin/env nextflow


  ch_ont = Channel
        .fromPath(params.ont)
        .view{"Input ONT: $it"}

  ch_illumina = Channel
        .fromFilePairs(params.illumina)
        .view{"Input Illumina: $it"}

  ch_illumina_samplelist = Channel.empty()
  ch_ont_fastq = Channel.empty()
  ch_ont_index = Channel.empty()
  ch_ont_ids = Channel.empty()
  ch_aligned = Channel.empty()
  ch_pcs = Channel.empty()
  ch_pcs_all = Channel.empty()
  ch_binned = Channel.empty()
  ch_assigned = Channel.empty()
  output_ch = Channel.empty()

process getIlluminaSampleList{
    label 'nf_01_ilslst'
    cpus params.resources.standard1.cpus
    memory params.resources.standard1.mem
    errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
    maxRetries 5
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

process makeFastaFromFastq {
    label 'nf_02_gtfa'
    cpus params.resources.standard2.cpus
    memory params.resources.standard2.mem
    errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
    maxRetries 5
    publishDir "$results_dir/2_makeFastaFromFastq", mode: 'symlink'
    input:
    path ont_file
    
    output:
    path('*.fasta.gz')

    shell:
    '''
    #Get fasta from nanopore fastq
    outfile=$(basename -s .fastq.gz !{ont_file} | sed "s/_/-/g")
    seqtk seq -a !{ont_file} > $outfile'.fasta'
    gzip $outfile'.fasta'

    '''

}

process getFastaIDs {
    label 'nf_03_fids'
    cpus params.resources.standard1.cpus
    memory params.resources.standard1.mem
    publishDir "$results_dir/3_getFastaIDs", mode: 'symlink'
    errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
    maxRetries 5
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

process IndexReferenceBwa {
    label 'nf_04_idx'
    cpus params.resources.standard1.cpus
    memory params.resources.standard1.mem
    errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
    maxRetries 5
    publishDir "$results_dir/4_ont_index_bwa", mode: 'symlink'

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
  label 'nf_05_algn'
  cpus params.resources.alignment.cpus
  memory params.resources.alignment.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  publishDir "$results_dir/5_alignIllumina", mode: 'symlink'
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
  label 'nf_06_smflt'
  cpus params.resources.samtoolsfilter.cpus
  memory params.resources.samtoolsfilter.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  publishDir "$results_dir/6_filterIlluminaAlignment", mode: 'symlink'
  input:
    path bam
    val mapq
    val include_flag_f
    val exclude_flag_F

  output:
    path('*.bam')

  shell:

  '''
      outfile=$(basename -s _raw.bam !{bam})
      filtbam=$outfile'.bam'

      #Filter and sort sam alignment file then output as bam
      samtools view -@ !{params.resources.samtoolsfilter.cpus} -bu -F !{exclude_flag_F} -f !{include_flag_f} -q !{mapq} !{bam} | samtools sort -o $filtbam
      samtools index $filtbam
  '''

}
//For now skipped, samtools depth performed in calculatePercentCovered and python script reads from stdin
process getCoverage{
  label 'nf_07_dpth'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  publishDir "$results_dir/7_getCoverage", mode: 'copy'
  input:
    path bam

  output:
    path('*.depth.gz')

  shell:
  '''
      outfile=$(basename -s .bam !{bam})
      samtoolsdepth=$outfile'.depth'

      #Get coverage at each position on each nanopore read
      samtools depth -a !{bam} > $samtoolsdepth
      gzip $samtoolsdepth
  '''
}

process calculatePercentCovered{
  label 'nf_08_pccov'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  publishDir "$results_dir/8_calculatePercentCovered", mode: 'symlink'
  input:
    path bam

  output:
    path '*.pc'

  shell:
  '''
      outfile=$(basename -s .bam !{bam})
      pcfile=$outfile'.pc'

      #Calculate percent coverage for each nanopore read
      samtools depth -a !{bam} | python !{params.scriptsdir}percent_coverage_from_depth_file.py -i stdin -f $pcfile > $pcfile

  '''

}

process concatenatePC {
   label 'nf_09_catpc'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  publishDir "$results_dir/9_concatenatePC", mode: 'copy'
  input:
    tuple(val(ont_id), val(pcs))

  output:
    path "${ont_id}_all.pc"

  shell:
  '''
  for pc in !{pcs}
  do
    echo $pc >> TMP.TXT
    pc=$(tr -d ,[] <<<$pc )
    cat $pc >> !{ont_id}_all.pc
  done
  '''
}

process binOntReadsToSpecies {
  label 'nf_10_br2sp'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  publishDir "$results_dir/10_binOntReadsToSpecies", mode: 'copy'
  input:
    tuple(path(all_pcs), path(ont_ids))
  output:
    path "*.binned"
  shell:
  '''
  #Bin each nanopore read to reference which has the highest pc

  outfile=$(basename -s .pc !{all_pcs})
  python !{params.scriptsdir}minion_read_bin_from_perc_cov.py !{all_pcs} !{ont_ids}  > $outfile'.binned'

  '''
}

process countReadsPerReference{
  label 'nf_11_crpr'
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 119..120 ? 'retry' : 'terminate' }
  maxRetries 5
  storeDir "$results_dir/11_countReadsPerReference"
  input:
    path binned_reads
    path illumina_ids
    val min_perc
    val max_perc
    val threshold_pct
  output:
    file "*_bin_counts.tsv"
shell:
  '''
  #Count number of reads binned to each reference and calculate percentages
  echo !{binned_reads} >> TMP.TXT
  outfile=$(basename -s .binned !{binned_reads})
  python !{params.scriptsdir}minion_read_counts_and_pcts.py !{binned_reads} !{illumina_ids} $outfile'_bin_counts.tsv'  !{min_perc} !{max_perc} !{threshold_pct}
  '''
}

workflow {

  makeFastaFromFastq(ch_ont)
  ch_ont_fastq = makeFastaFromFastq.out
  ch_ont_fastq.view{ "Fasta created from fastq: $it" }

  IndexReferenceBwa(ch_ont_fastq)
  ch_ont_index = IndexReferenceBwa.out
  ch_ont_index.view{ "BWA index created from fasta: $it" }

  getFastaIDs(ch_ont_fastq)
  ch_ont_ids = getFastaIDs.out
  ch_ont_ids.view{ "Fasta ID list created from fasta: $it" }

  ch_illumina_samples = ch_illumina.map{x ->
      def sname = x.get(0)
      return sname}
    .collect()
    .view()
  getIlluminaSampleList(ch_illumina_samples)
  ch_illumina_samplelist=getIlluminaSampleList.out

  ch_ont_index=ch_ont_index.combine(ch_illumina)
    .view{ "ch_ont_index combined: $it" }
  
  
  alignIllumina(ch_ont_index) 
  ch_aligned = alignIllumina.out.view{ "Alignment result: $it" }
  filterIlluminaAlignment(
    ch_aligned,
    params.filterIlluminaAlignment.mapq,
    params.filterIlluminaAlignment.include_flag_f,
    params.filterIlluminaAlignment.exclude_flag_F
  ) | 
  view{ "filterIlluminaAlignmentignment result: $it" } |
  //getCoverage | view{ "getCoverage result: $it" }  |
  calculatePercentCovered | view { "calculatePercentCovered result: $it" }
  
  //Percent of each read that is covered by each species
  ch_pcs = calculatePercentCovered.out.map{ file ->
      def key = file.name.toString().tokenize('_').get(0)
      return tuple(key, file)
    }
    .groupTuple()
    .view{ "calculatePercentCovered grouped by nanopore sample: $it" }

  //Merge all pcs from one nanopore sample into one file
  concatenatePC(ch_pcs)
  
  //Combine each nanopopre sample pc with its read ids
  ch_pcs_all = concatenatePC.out.view{ "concatenatePC result: $it" }
           .phase(ch_ont_ids) { it -> 
              it.name.toString().replaceFirst(/_all.pc/, ".ids")
            }.view{"Phase .pc and .fasta.ids: $it"}
  binOntReadsToSpecies(ch_pcs_all)
  ch_binned = binOntReadsToSpecies.out.view{"binOntReadsToSpecies output: $it"}
  countReadsPerReference(ch_binned, 
                        ch_illumina_samplelist, 
                        params.countReadsPerReference.min_perc, 
                        params.countReadsPerReference.max_perc,
                        params.countReadsPerReference.threshold_pct)
  ch_assigned = countReadsPerReference.out.view{"countReadsPerReference output: $it"}

}