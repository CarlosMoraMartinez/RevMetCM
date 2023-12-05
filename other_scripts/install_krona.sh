#previous:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

#create krona_env with krona
conda create -n krona_env krona  -y
conda activate krona_env

# delete a symbolic link that is not correct
rm -rf /home/carmoma/anaconda3//envs/krona_env/opt/krona/taxonomy

#  create a directory in  home where the krona database will live
mkdir -p ~/krona/taxonomy

# now we make a symbolic link to that directory
ln -s ~/krona/taxonomy /home/carmoma/anaconda3/envs/krona_env/opt/krona/taxonomy

# the krona downloader fails when downloading the latest NCBI db
# ktUpdateTaxonomy.sh ~/krona/taxonomy
# taxdump.tar.gz needs to be downloaded manually from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
# taxdump.tar.gz needs to be placed in ~/krona/taxonomy
# then one needs to run  ktUpdateTaxonomy.sh --only-build
ktUpdateTaxonomy.sh --only-build

cd ~/projects/pollen/results/problem_samples_1/kraken_taxons
mkdir krona_charts
for i in $(ls *txt.gz); do echo $i; ktImportTaxonomy -m 3 -t 5 <(zcat $i) -o krona_charts/${i%.txt.gz}.html; done
ktImportTaxonomy -m 3 -t 5 <(zcat *.txt.gz) -o krona_charts/all_mock.html

cd ~/projects/pollen/results/mock_kraken1/kraken_tx_ont
cd ~/projects/pollen/results/mock_kraken1/kraken_tx_illumina


##PUCCINIA STRIIFORMIS: patógeno del trigo encontrado en polen mock samples.
