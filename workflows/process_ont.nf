include { filterOntReads } from '../modules/filter_ont_reads'
include { makeFastaFromFastq } from '../modules/make_fasta_from_fastq'
include { getFastaIDs } from '../modules/get_fasta_ids'
include { indexReference } from '../modules/index_reference'
include { getLowComplexityRegions } from '../modules/get_low_complexity_regions'
include { taxonFilterOnt } from '../modules/taxon_filter_ont'
include { organelleGenomeIndex } from '../modules/organelle_genome_index'
include { organelleFilterOnt } from '../modules/organelle_filter_ont'


workflow ONT2FASTA {
  take: ch_ont
  main:
    if(params.filterOntReads.do_filter){
      filterOntReads(ch_ont, params.filterOntReads.min_length, params.filterOntReads.min_mean_q)
      ch_ont = filterOntReads.out
      //ch_ont.view{ "Fastq filtered filtlong: $it" }
    }

    ch_taxons_ont = Channel.from([])
    if(params.taxonFilterOnt.do_filter){
      taxonFilterOnt(ch_ont, params.taxonFilterOnt.db)
      ch_ont = taxonFilterOnt.out.map{it -> it[1]}
        //.view{ "ONT fastq filtered by taxon: $it" }
      ch_taxons_ont = taxonFilterOnt.out.map{it -> it[0]}
        //.view{ "ONT taxon composition: $it" }
    }
    if(params.organelleFilter.do_filter){
      if(params.organelleFilter.from_fasta){
        organelleGenomeIndex(params.organelleFilter.db)
        ch_organelle_index = organelleGenomeIndex.out
      }else{
        ch_organelle_index = params.organelleFilter.index
      }
      //ch_organelle_index.combine(ch_organelle_index, ch_ont)
      //  .view{ "ONT organelle filter input: $it" }
      organelleFilterOnt(ch_organelle_index, ch_ont)
      ch_ont = organelleFilterOnt.out
      .view{ "ONT organelle filter output: $it" }
    }

    makeFastaFromFastq(ch_ont)
    ch_ont_fasta = makeFastaFromFastq.out
    //ch_ont_fasta.view{ "Fasta created from fastq: $it" }

    indexReference(ch_ont_fasta)
    ch_ont_index = indexReference.out
    //ch_ont_index.view{ "BWA index created from fasta: $it" }
    getFastaIDs(ch_ont_fasta)
    ch_ont_ids = getFastaIDs.out 
    //ch_ont_ids.view{ "Fasta ID list created from fasta: $it" }

    ch_ont_dust = Channel.from([])
    if(params.filterLowComplexityRegions){
      getLowComplexityRegions(ch_ont_fasta,
        params.getLowComplexityRegions.level, 
        params.getLowComplexityRegions.window)
      ch_ont_dust = getLowComplexityRegions.out
      //ch_ont_dust.view{ "BED with low complexity regions: $it" }
    }
    
  emit:
    ch_ont_index
    ch_ont_ids
    ch_taxons_ont
    ch_ont_dust
    ch_ont
}