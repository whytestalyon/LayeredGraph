#!/usr/bin/env bash
# WARNING this script is not designed to run on its own
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

#pubtator parser
echo "Downloading HPO disease to phenotype associations..."
curl --retry 3 -L http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt -o ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt
echo "Downloading MedGen data..."
curl --retry 3 -L ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/MedGen_HPO_OMIM_Mapping.txt.gz -o MedGen_HPO_OMIM_Mapping.txt.gz
echo "Downloading non-alt loci gene info from the HGNC..."
curl --retry 3 -L ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt -o non_alt_loci_set.txt
echo "Downloading the Disease Ontology..."
curl --retry 3 -L http://purl.obolibrary.org/obo/doid.obo -o doid.obo
echo "Downloading and unpacking OrphaNet diseseas ontology..."
curl --retry 3 -L  http://www.orphadata.org/data/ORDO/ordo_orphanet.owl.zip -o ordo_orphanet.owl.zip
echo "Downloading PubTator bioconcepts..."
curl --retry 3 -L ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator/bioconcepts2pubtator.gz -o bioconcepts2pubtator.gz

#ppi data files
echo "Downloading PPI graph files..."
curl --retry 3 -g -L "http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes&force=no&preview=true&format=tab" -o uniprot-all.tab.gz
gunzip uniprot-all.tab.gz
curl --retry 3 -L "http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz" -o ConsensusPathDB_human_PPI.gz
gunzip ConsensusPathDB_human_PPI.gz

cd ..
#we need the HPO to be included
cp ./HPO_data_files/hp.obo ./HPO_graph_data/hp.obo

# run the pubtator parser and generate the gene to phenotype correlations from it
echo
echo "Parsing PubTator data:"
python3 -u ./PubTator/PubTatorParser.py
echo

#we need to generate the HPO graph data
echo
echo "Building primary HPO graph:"
python3 -u ./LayeredGraphAPI/PyxisMapBuilder.py

#we need to generate the PPI graph data
echo
echo "Building primary PPI graph:"
python3 -u ./LayeredGraphAPI/ProtBuilder.py
