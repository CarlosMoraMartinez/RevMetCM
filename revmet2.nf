#!/usr/bin/env nextflow

params.ont = "/home/carmoma/projects/pollen/revmet_paper_data/ena_files_nanopore/ERR3132323/mock_mixes_nanopore_reads/mock_mixes/barcode*.fastq"
params.illumina = "/home/carmoma/projects/pollen/revmet_paper_data/ena_files_illumina/*_R{1,2}*.fastq.gz"

  ch_ont = Channel
        .fromPath(params.ont)
        .view{"Input ONT: $it"}

  ch_illumina = Channel
        .fromFilePairs(params.illumina)
        .view{"Input Illumina: $it"}

  
  ch_illumina_samplelist = Channel.empty()
  ch_align = Channel.empty()
  ch_ont_ids = Channel.empty()
  ch_pcs = Channel.empty()
  ch_pcs_all = Channel.empty()
  ch_binned = Channel.empty()


process getIlluminaSampleList{
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

    input:
    tuple(path(ont_file), val(illumina_id), val(illumina_reads))

    
    output:
    path '*.pc'
    
    shell:
    '''
    outfile=$(basename -s .fasta !{ont_file})_!{illumina_id}
    echo 'ILL ID: !{illumina_id}, FICHEROS: !{illumina_reads[0]} - !{illumina_reads[1]} -- !{ont_file} ' >> $outfile'.pc'
    '''

}

process concatenatePC {
  input:
  tuple(val(ont_id), val(pcs))

  output:
  path "all_${ont_id}.pc"

  shell:
  '''
  for pc in !{pcs}
  do
    pc=$(tr -d ,[] <<<$pc )
    cat $pc >> all_!{ont_id}.pc
  done
  '''
}

process binOntReadsToSpecies {
  input:
    path all_pcs
    path ont_ids
  output:
   path "*.binned"
  shell:
  '''
  #Bin each nanopore read to reference which has the highest pc
  outfile=$(basename -s .pc !{all_pcs})
  python /home/carmoma/projects/pollen/myrevmet/mock.py !{all_pcs} !{ont_ids}  > $outfile'.binned'

  '''
}

workflow {

  makeFastaFromFastq(ch_ont)
  ch_align = makeFastaFromFastq.out
  ch_align.view{ "Fasta created from fastq: $it" }

  getFastaIDs(ch_align)
  ch_ont_ids = getFastaIDs.out
  ch_ont_ids.view{ "Fasta ID list created from fasta: $it" }

  ch_illumina_samples = ch_illumina.map{x ->
      def sname = x.get(0)
      return sname}
    .collect()
    .view()
  getIlluminaSampleList(ch_illumina_samples)
  ch_illumina_samplelist=getIlluminaSampleList.out

  ch_align=ch_align.combine(ch_illumina)
    .view{ "ch_align combined: $it" }
  alignIllumina(ch_align)
  ch_pcs = alignIllumina.out.map{ file ->
      def key = file.name.toString().tokenize('_').get(0)
      return tuple(key, file)
    }
    .groupTuple()
    .view{ "Alignment result: $it" }

  concatenatePC(ch_pcs)
  ch_pcs_all = concatenatePC.out.view()

  binOntReadsToSpecies(ch_pcs_all, ch_ont_ids)
  ch_binned = binOntReadsToSpecies.out.view()
}