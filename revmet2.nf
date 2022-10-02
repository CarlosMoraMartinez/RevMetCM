#!/usr/bin/env nextflow


  ch_ont = Channel
        .fromPath(params.ont)
        .view{"Input ONT: $it"}

  ch_illumina = Channel
        .fromFilePairs(params.illumina)
        .view{"Input Illumina: $it"}
  
  //results_dir = Channel.value(params.results_dir).view{"Output dir: $it"}

  ch_illumina_samplelist = Channel.empty()
  ch_ont_fastq = Channel.empty()
  ch_ont_ids = Channel.empty()
  ch_aligned = Channel.empty()
  ch_pcs = Channel.empty()
  ch_pcs_all = Channel.empty()
  ch_binned = Channel.empty()
  ch_assigned = Channel.empty()
  output_ch = Channel.empty()

process getIlluminaSampleList{
    publishDir "$results_dir/1_getIlluminaSampleList", mode: 'symlink'
    input:
    val ill_names
    
    output:
    path 'skim_ref.ids'

    shell:
    '''
    echo !{ill_names} | sed 's/, /\\n/g' | tr -d [] > skim_ref.ids

    '''

}
process makeFastaFromFastq {
    publishDir "$results_dir/2_makeFastaFromFastq", mode: 'symlink'
    input:
    path ont_file
    
    output:
    path('*.fasta')

    shell:
    '''
    #Get fasta from nanopore fastq
    outfile=$(basename -s .fastq !{ont_file})
    seqtk seq -a !{ont_file} > $outfile'.fasta'

    '''

}

process getFastaIDs {
    publishDir "$results_dir/3_getFastaIDs", mode: 'symlink'

    input:
    path fasta_file
    
    output:
     path('*.ids')


    shell:
    '''
    #Get list of nanopore read ids from nanopore fasta
    outfile=$(basename -s .fastq !{fasta_file})
    grep ">" !{fasta_file}  | sed 's/>//g' | cut -f 1 -d\\  > $outfile'.ids'
    '''

}

process alignIllumina {
  publishDir "$results_dir/4_alignIllumina", mode: 'symlink'
  input:
    tuple(path(ont_file), val(illumina_id), val(illumina_reads))

    
  output:
    path('*_raw.sam')
    
  shell:
  if( params.alignIllumina.program == 'minimap2' )
      '''
        outfile=$(basename -s .fasta !{ont_file})_!{illumina_id}
        rawsam=$outfile'_raw.sam'
  
        minimap2 -ax sr !{ont_file} !{illumina_reads[0]}  !{illumina_reads[1]} > $rawsam
      '''
  else if( params.alignIllumina.program == 'bwa' )
      '''
        outfile=$(basename -s .fasta !{ont_file})_!{illumina_id}
        rawsam=$outfile'_raw.sam'
        
        bwa index !{ont_file}
        bwa mem !{ont_file} !{illumina_reads[0]}  !{illumina_reads[1]} > $rawsam
      '''


}

process filterIlluminaAlignment{
  publishDir "$results_dir/5_filterIlluminaAlignment", mode: 'symlink'
  input:
    path sam
    val mapq
    val include_flag_f
    val exclude_flag_F

  output:
    path('*.bam')

  shell:

  '''
      outfile=$(basename -s _raw.sam !{sam})
      filtbam=$outfile'.bam'

      #Filter and sort sam alignment file then output as bam
      samtools view -bu -F !{exclude_flag_F} -f !{include_flag_f} -q !{mapq} !{sam} | samtools sort -o $filtbam
      samtools index -c $filtbam
  '''

}

process getCoverage{
  publishDir "$results_dir/6_getCoverage", mode: 'symlink'
  input:
    path bam

  output:
    path('*.depth')

  shell:
  '''
      outfile=$(basename -s .bam !{bam})
      samtoolsdepth=$outfile'.depth'

      #Get coverage at each position on each nanopore read
      samtools depth -a !{bam} > $samtoolsdepth
  '''
}

process calculatePercentCovered{
  publishDir "$results_dir/7_calculatePercentCovered", mode: 'symlink'
  input:
    path depthfile

  output:
    path '*.pc'

  shell:
  '''
      outfile=$(basename -s .sam !{depthfile})
      pcfile=$outfile'.pc'

      #Calculate percent coverage for each nanopore read
      python /home/carmoma/projects/pollen/myrevmet/scripts/percent_coverage_from_depth_file.py !{depthfile} > $pcfile

  '''

}

process concatenatePC {
  publishDir "$results_dir/8_concatenatePC", mode: 'symlink'
  input:
    tuple(val(ont_id), val(pcs))

  output:
    path "${ont_id}_all.pc"

  shell:
  '''
  for pc in !{pcs}
  do
    pc=$(tr -d ,[] <<<$pc )
    cat $pc >> !{ont_id}_all.pc
  done
  '''
}

process binOntReadsToSpecies {
  publishDir "$results_dir/9_binOntReadsToSpecies", mode: 'symlink'
  input:
    tuple(path(all_pcs), path(ont_ids))
  output:
    path "*.binned"
  shell:
  '''
  #Bin each nanopore read to reference which has the highest pc

  outfile=$(basename -s .pc !{all_pcs})
  python /home/carmoma/projects/pollen/myrevmet/scripts/minion_read_bin_from_perc_cov.py !{all_pcs} !{ont_ids}  > $outfile'.binned'

  '''
}

process countReadsPerReference{
  publishDir "$results_dir/10_countReadsPerReference", mode: 'symlink'
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

  outfile=$(basename -s .binned !{binned_reads})
  python /home/carmoma/projects/pollen/myrevmet/scripts/minion_read_counts_and_pcts.py !{binned_reads} !{illumina_ids} $outfile'_bin_counts.tsv'  !{min_perc} !{max_perc} !{threshold_pct}
  '''
}

workflow {

  makeFastaFromFastq(ch_ont)
  ch_ont_fastq = makeFastaFromFastq.out
  ch_ont_fastq.view{ "Fasta created from fastq: $it" }

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

  ch_ont_fastq=ch_ont_fastq.combine(ch_illumina)
    .view{ "ch_ont_fastq combined: $it" }
  
  
  alignIllumina(ch_ont_fastq) 
  ch_aligned = alignIllumina.out.view{ "Alignment result: $it" }
  filterIlluminaAlignment(
    ch_aligned,
    params.filterIlluminaAlignment.mapq,
    params.filterIlluminaAlignment.include_flag_f,
    params.filterIlluminaAlignment.exclude_flag_F
  ) | 
  view{ "filterIlluminaAlignmentignment result: $it" } |
  getCoverage | view{ "getCoverage result: $it" }  |
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
  ch_pcs_all = concatenatePC.out.view()
           .phase(ch_ont_ids) { it -> 
              it.name.toString().replaceFirst(/_all.pc/, ".fasta.ids")
            } .view{"Phase .pc and .fasta.ids: $it"}
  binOntReadsToSpecies(ch_pcs_all)
  ch_binned = binOntReadsToSpecies.out.view{"binOntReadsToSpecies output: $it"}
  countReadsPerReference(ch_binned, 
                        ch_illumina_samplelist, 
                        params.countReadsPerReference.min_perc, 
                        params.countReadsPerReference.max_perc,
                        params.countReadsPerReference.threshold_pct)
  ch_assigned = countReadsPerReference.out.view{"countReadsPerReference output: $it"}

}