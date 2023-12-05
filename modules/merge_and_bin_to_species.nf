process mergeAndBin2species {
  label 'rvm05_bin2sp'
  conda params.mergeAndBin2species.conda
  cpus params.resources.mergeandbin2species.cpus
  memory params.resources.mergeandbin2species.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 3
  storeDir "$results_dir/rvm05_MergeAndBin"

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
