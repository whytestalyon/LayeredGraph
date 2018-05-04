#!/usr/bin/env bash
# WARNING this script is not designed to run on its own
#pull all files
mkdir HPO_graph_data
mkdir HPO_data_files

cd HPO_data_files

if [[ "download" == "$1" || ! -e hp.obo ]]; then
    echo "Getting latest HPO graph structure..."
    curl --retry 3 -L http://purl.obolibrary.org/obo/hp.obo -o hp.obo
fi

if [[ "download" == "$1" || ! -e ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt ]]; then
    echo "Getting latest HPO to gene relationships..."
    curl --retry 3 -L http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt -o ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
fi


if [[ "download" == "$1" || ! -e hgnc_complete_set.txt ]]; then
    echo "Retrieving latest HGNC complete set..."
    curl --retry 3 -L ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt -o hgnc_complete_set.txt
fi

if [[ "download" == "$1" || ! -e ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt ]]; then
    echo "Downloading HPO disease to phenotype associations..."
    curl --retry 3 -L http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt -o ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt
fi


if [[ "download" == "$1" || ! -e MedGen_HPO_OMIM_Mapping.txt.gz ]]; then
    echo "Downloading MedGen data..."
    curl --retry 3 -L ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/MedGen_HPO_OMIM_Mapping.txt.gz -o MedGen_HPO_OMIM_Mapping.txt.gz
fi


if [[ "download" == "$1" || ! -e non_alt_loci_set.txt ]]; then
    echo "Downloading non-alt loci gene info from the HGNC..."
    curl --retry 3 -L ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt -o non_alt_loci_set.txt
fi


if [[ "download" == "$1" || ! -e doid.obo ]]; then
    echo "Downloading the Disease Ontology..."
    curl --retry 3 -L http://purl.obolibrary.org/obo/doid.obo -o doid.obo
fi


if [[ "download" == "$1" || ! -e ordo_orphanet.owl.zip ]]; then
    echo "Downloading and unpacking OrphaNet diseseas ontology..."
    curl --retry 3 -L  http://www.orphadata.org/data/ORDO/ordo_orphanet.owl.zip -o ordo_orphanet.owl.zip
fi


if [[ "download" == "$1" || ! -e gene2pubtator.gz ]]; then
    echo "Downloading PubTator gene file..."
    curl --retry 3 -L ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator/gene2pubtator.gz -o gene2pubtator.gz
fi


if [[ "download" == "$1" || ! -e disease2pubtator.gz ]]; then
    echo "Downloading PubTator disease file..."
    curl --retry 3 -L ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator/disease2pubtator.gz -o disease2pubtator.gz
fi

#ppi data files
if [[ "download" == "$1" || ! -e uniprot-all.tab ]]; then
    echo "Downloading PPI graph files..."
    curl --retry 3 -g -L "http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes&force=no&preview=true&format=tab" -o uniprot-all.tab.gz
    gunzip uniprot-all.tab.gz
fi


if [[ "download" == "$1" || ! -e ConsensusPathDB_human_PPI ]]; then
    curl --retry 3 -L "http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz" -o ConsensusPathDB_human_PPI.gz
    gunzip ConsensusPathDB_human_PPI.gz
fi

cd ..

#we need the HPO to be included
cp ./HPO_data_files/hp.obo ./HPO_graph_data/hp.obo

# run the pubtator parser and generate the gene to phenotype correlations from it
echo
echo "Parsing PubTator data:"
python3 -u ./PubTator/PubTatorParser.py
echo

# need the pubtator output included
cp ./HPO_data_files/gene2phenotype.json ./HPO_graph_data/gene2phenotype.json
cp ./HPO_data_files/gene2phenotypeIndex.json ./HPO_graph_data/gene2phenotypeIndex.json
cp ./HPO_data_files/non_alt_loci_set.txt ./HPO_graph_data/non_alt_loci_set.txt

#we need to generate the HPO graph data
echo
echo "Building primary HPO graph:"
python3 -u ./LayeredGraphAPI/PyxisMapBuilder.py

#we need to generate the PPI graph data
echo
echo "Building primary PPI graph:"
python3 -u ./LayeredGraphAPI/ProtBuilder.py

#build the metadata file
echo
echo "Building the metadata file:"
python3 -u ./buildMetadata.py
