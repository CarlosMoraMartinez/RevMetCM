include { alignIllumina } from '../modules/align_illumina'
include { filterIlluminaAlignment } from '../modules/filter_illumina_alignment'
include { filterDustRegions } from '../modules/filter_dust_regions'
include { getCoverage } from '../modules/get_cogerage'
include { mergeAndBin2species } from '../modules/merge_and_bin_to_species'


workflow REVMETCORE {
  take: 
    ch_ont_index
    ch_ont_dust
    ch_illumina_processed

  main:
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
      ch_dustbed = ch_ont_dust
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
    ch_merged = mergeAndBin2species.out

  emit:
    ch_merged
}