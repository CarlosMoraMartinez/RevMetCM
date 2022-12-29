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

workflow illuminafastq {
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
  emit:
    ch_illumina_samplelist
}