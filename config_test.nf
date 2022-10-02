
basedir = "/home/carmoma/projects/pollen/revmet_paper_data/"
params.ont = "$basedir/ena_files_nanopore/ERR3132323/mock_mixes_nanopore_reads/mock_mixes/barcode*.fastq"
params.illumina = "$basedir/ena_files_illumina/cortados/*_R{1,2}*.fastq.gz"

env.results_dir = "/home/carmoma/projects/pollen/myrevmet/resultados_test2"

params.mapping.program = "bwa"

params.filterIlluminaAlignment.include = 0
//read unmapped, not primary alignment, supplementary alignment
params.filterIlluminaAlignment.exclude = 2308
params.filterIlluminaAlignment.mapq = 0

params.min_perc = 1.0
params.max_perc = 99.9

dag {
    enabled = true
    file = 'pipeline_dag.html'
}