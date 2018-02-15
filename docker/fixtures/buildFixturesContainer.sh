#!/usr/bin/env bash
#echo "Preparing to download graph data, please enter Morgan username:"
#read USERNAME
#echo "Downloading graph data from Morgan, scp will prompt for password..."
#scp -rp ${USERNAME}@login.morgan.haib.org:/gpfs/gpfs1/home/jholt/HPO_graph_data .

#pull all files
mkdir HPO_graph_data
cd HPO_graph_data
echo "Getting latest HPO graph structure..."
curl -L http://purl.obolibrary.org/obo/hp.obo -o hp.obo
echo "Getting latest HPO to gene relationships..."
curl -L http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt -o ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
cd ..

#we need to generate the HPO graph data


#we need to generate the PPI graph data


#docker build -t docker-registry.haib.org/sdi/layered-graph-fixtures .
#if [ $? -ne 0 ]; then
#    echo "Failed to build the fixtures container..."
#else
#    echo "Fixtures image built!"
#fi

#echo "Cleaning up..."
#rm -rf HPO_graph_data

