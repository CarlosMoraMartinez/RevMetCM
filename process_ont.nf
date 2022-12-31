#!/usr/bin/env nextflow

process filterOntReads {
  label 'ont01_flng'
  conda params.filterOntReads.conda
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont01_filterOntReads-l$min_length-q$min_mean_q", mode: 'copy'
  
  input:
  path ont_file
  val min_length
  val min_mean_q
  
  output:
  path('*trim.fastq.gz')

  shell:
  '''
  #Get fasta from nanopore fastq
  outfile=$(basename -s .fastq.gz !{ont_file} | sed "s/_/-/g")
  filtlong --min_length !{min_length} --min_mean_q !{min_mean_q} !{ont_file} | gzip > $outfile'.trim.fastq.gz'
  '''

}

process makeFastaFromFastq {
  label 'ont02_gtfa'
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont02_makeFastaFromFastq", mode: 'symlink'
  
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
  label 'ont03_fids'
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  publishDir "$results_dir/ont03_getFastaIDs", mode: 'symlink'
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

process indexReference {
  label 'ont04_idx'
  cpus params.resources.index.cpus
  memory params.resources.index.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont04_index", mode: 'symlink'

  input:
  path fasta_file
  
  output:
  tuple path("${fasta_file}"), path("${fasta_file}*")


  shell:
  if(params.alignIllumina.program == "minimap2")
    '''
    minimap2 -x sr -d !{fasta_file}'.mmi' !{fasta_file}
    '''
  else if(params.alignIllumina.program == "bwa")
    '''
    #Index nanopore fasta reference prior to alignment
    #Note: this will happen even if minimap2 is used instead of bwa. Fix in the future.
    bwa index !{fasta_file}
  '''
}


process getLowComplexityRegions {
  label 'ont05_dust'
  conda params.getLowComplexityRegions.conda
  cpus params.resources.standard2.cpus
  memory params.resources.standard2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont05_dust", mode: 'symlink'

  input:
  path fasta_file
  val window
  val level
  
  output:
  path('*.bed')

  shell:
  '''
  outfile=$(basename -s .fasta.gz !{fasta_file})
  dustmasker -in <(zcat !{fasta_file}) -outfmt acclist -level !{level} -window !{window} > $outfile'.intlist'
  paste <(awk '{print $1}' $outfile'.intlist') <(cut -f 2,3 $outfile'.intlist') | sed "s/>//" > $outfile'.bed'
  rm $outfile'.intlist'
  '''
}

process taxonFilterOnt {
  label 'ont06_taxonFiltOnt'
  conda params.taxonFilterOnt.conda
  cpus params.resources.taxonFilterOnt.cpus
  memory params.resources.taxonFilterOnt.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ont06_taxonFiltOnt", mode: 'symlink'

  input:
  path fastq
  path db
  
  output:
  tuple(path("*taxons.txt.gz"), path("*.fastq.gz"))

  shell:
  '''
  outfile=$(basename -s .fastq.gz !{fastq})
  
  kraken2 --db !{db} --threads !{params.resources.taxonFilterOnt.cpus} --unclassified-out $outfile'.tx.fastq' --gzip-compressed !{fastq} | gzip > $outfile'_taxons.txt.gz' 
  gzip $outfile'.tx.fastq' #TO DO: include pigz in the conda environment

  '''
}

workflow ont2fasta {
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

    makeFastaFromFastq(ch_ont)
    ch_ont_fasta = makeFastaFromFastq.out
    ch_ont_fasta.view{ "Fasta created from fastq: $it" }

    indexReference(ch_ont_fasta)
    ch_ont_index = indexReference.out
    ch_ont_index.view{ "BWA index created from fasta: $it" }
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