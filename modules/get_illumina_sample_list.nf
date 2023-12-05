process getIlluminaSampleList{
  label 'ilm01_ilslst'
  conda params.getIlluminaSampleList.conda
  cpus params.resources.standard1.cpus
  memory params.resources.standard1.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/ilm01_getIlluminaSampleList"
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