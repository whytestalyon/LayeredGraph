#!/usr/bin/env bash
#pull all files
mkdir HPO_graph_data
mkdir HPO_data_files
cd HPO_data_files
echo "Getting latest HPO graph structure..."
curl --retry 3 -L http://purl.obolibrary.org/obo/hp.obo -o hp.obo
echo "Getting latest HPO to gene relationships..."
curl --retry 3 -L http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt -o ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
echo "Retrieving latest HGNC complete set..."
curl --retry 3 -L ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt -o hgnc_complete_set.txt

#temporary until I get source files from Brandon
echo "Downloading temporary PubTator files, enter Morgan username:"
read USERNAME
echo "Downloading data, scp will prompt for password..."
scp -rp ${USERNAME}@login.morgan.haib.org:/gpfs/gpfs1/home/jholt/temp/gene2phenotype.json .

#ppi data files
curl --retry 3 -g -L "http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes&force=no&preview=true&format=tab" -o uniprot-all.tab.gz
gunzip uniprot-all.tab.gz
curl --retry 3 -L "http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz" -o ConsensusPathDB_human_PPI.gz
gunzip ConsensusPathDB_human_PPI.gz

cd ..

#we need the HPO to be included
cp ./HPO_data_files/hp.obo ./HPO_graph_data/hp.obo

#we need to generate the HPO graph data
echo -e "\nBuilding primary HPO graph:"
python3 ../../LayeredGraphAPI/PyxisMapBuilder.py

#we need to generate the PPI graph data
echo -e "\nBuilding primary PPI graph:"
python3 ../../LayeredGraphAPI/ProtBuilder.py

ls -l HPO_graph_data

docker build -t docker-registry.haib.org/sdi/layered-graph-fixtures .
if [ $? -ne 0 ]; then
    echo "Failed to build the fixtures container..."
else
    echo "Fixtures image built!"
fi

echo "Cleaning up..."
#uncomment for local dev
cp -r ./HPO_graph_data ../../HPO_graph_data
rm -rf HPO_data_files
rm -rf HPO_graph_data
