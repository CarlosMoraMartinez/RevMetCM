
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




