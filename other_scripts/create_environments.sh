#previous:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict


conda create -y --name fastqc-env python=3.10
conda activate fastqc-env
conda install -y -c bioconda fastqc

conda create -y --name trimmomatic-env python=3.10
conda activate trimmomatic-env
conda install -y -c bioconda trimmomatic

conda create -y --name bowtie2-env python=3.10
conda activate bowtie2-env
conda install -y -c bioconda bowtie2 samtools

conda create -y --name bwa-env python=3.10
conda activate bwa-env
conda install -y -c bioconda bwa bwa-mem2 samtools

#conda create -y --name kraken-env python=3.10
#conda activate kraken-env
#conda install -y -c bioconda kraken 
#conda install -y -c bioconda krakentools 

conda create -y --name kraken2-env python=3.10
conda activate kraken2-env
conda install -y -c bioconda kraken2
conda install -y -c bioconda pigz 

conda create -y --name bracken-env python=3.10
conda activate bracken-env
conda install -y -c bioconda bracken

conda create -y --name krona-env python=3.10
conda activate krona-env
conda install -y -c bioconda krona

conda create -y --name multiqc-env python=3.10
conda activate multiqc-env
#conda install -y -c bioconda multiqc
pip install multiqc

conda create -y --name seqtk-env python=3.10
conda activate seqtk-env
conda install -y -c bioconda seqtk pigz samtools

conda create -y --name seqkit-env python=3.10
conda activate seqkit-env
conda install -y -c bioconda seqkit pigz samtools

conda create -y --name humann-env python=3.10
conda activate humann-env
conda install -y -c bioconda humann

## Add Krakentools to kraken2-env
#get scripts from github
mkdir tmp
cd tmp
conda activate kraken2-env
conda install -y git
git clone https://github.com/jenniferlu717/KrakenTools

#get conda bin directory
krakendir=$(which kraken2)
condadir=$(dirname $krakendir)
echo $condadir

#Copy scripts to conda directory
cp -v KrakenTools/*.py $condadir 
chmod 777 $condadir/*.py

#Clean temporary directory
cd ../
rm -rf tmp

##########################
##Initialize Krona env
conda activate krona_env

# delete a symbolic link that is not correct
rm -rf /home/carmoma/anaconda3//envs/krona_env/opt/krona/taxonomy

#  create -y a directory in  home where the krona database will live
mkdir -p ~/krona/taxonomy

# now we make a symbolic link to that directory
ln -s ~/krona/taxonomy /home/carmoma/anaconda3/envs/krona_env/opt/krona/taxonomy

# the krona downloader fails when downloading the latest NCBI db
# ktUpdateTaxonomy.sh ~/krona/taxonomy
# taxdump.tar.gz needs to be downloaded manually from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
# taxdump.tar.gz needs to be placed in ~/krona/taxonomy
# then one needs to run  ktUpdateTaxonomy.sh --only-build
ktUpdateTaxonomy.sh --only-build

#Exampe commands for krona:
#for i in $(ls *txt.gz); do echo $i; ktImportTaxonomy -m 3 -t 5 <(zcat $i) -o krona_charts/${i%.txt.gz}.html; done
#ktImportTaxonomy -m 3 -t 5 <(zcat *.txt.gz) -o krona_charts/all_mock.html


## Nextflow environment (if needed to run the pipeline)

mkdir bin
cd bin
curl -s https://get.nextflow.io | bash
./nextflow run tutorial.nf

export PATH=$HOME/bin:$PATH
#And modify ~/.basrc file aadding this last ine


## R environment

conda create -y --name R-env python=3.10
conda activate R-env
conda install -c conda-forge r-base
conda install -c bioconda bioconductor-biocinstaller

## MASH environment
conda create -y --name mash-env python=3.10
conda activate mash-env
conda install -y -c bioconda mash