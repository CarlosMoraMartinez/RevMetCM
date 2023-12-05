
cd ~
miniconda3/bin/conda init bash
source .bashrc

conda activate /software/assembly/conda/kraken2.1.2/

#Standard:

kraken2-build --standard --threads 64 --db standard2

##Custom:

kraken2-build --download-taxonomy --db revmet

kraken2-build --download-library bacteria --db revmet
kraken2-build --download-library archaea --db revmet
kraken2-build --download-library fungi --db revmet
kraken2-build --download-library human --db revmet

kraken2-build --build --threads 64 --db revmet


#Errores al descargar librería:
#rsync_from_ncbi.pl: unexpected FTP path (new server?) for https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2

#Errores al ejecutar build custom:
#Creating sequence ID to taxonomy ID map (step 1)...
#No preliminary seqid/taxid mapping files found, aborting.

#/scratch/devel/talioto/denovo_assemblies/kraken_db/kraken2_comprehensive_DB_20181025/

#run
srun -c 128 --mem=300G --pty bash -li  

kraken2 --db kraken2_ArcBacPlaVirFunProUniVec_DB_20200714 --threads 120 --classified-out in.fastq --unclassified-out unclas.fastq --gzip-compressed ../dust_mock/2_filterOntReads-l1000-q7/FAR74611-1-NB01.trim.fastq.gz

kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

nohup kraken2 --db kraken2_comprehensive_DB_20181025 --threads 120 --classified-out in_compr#.fastq --unclassified-out unclas_compr#.fastq --gzip-compressed ../../easi_illumina2/AS0226.226.EASI_48.7161AF.HFYTVDRXY.2.210UDI-idt-UMI.1.fastq.gz  ../../easi_illumina2/AS0226.226.EASI_48.7161AF.HFYTVDRXY.2.210UDI-idt-UMI.2.fastq.gz

#A lot of kraken DB: https://benlangmead.github.io/aws-indexes/k2
# date 12/9/2022
wget https://benlangmead.github.io/aws-indexes/k2


