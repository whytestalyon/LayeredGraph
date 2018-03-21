#!/usr/bin/env bash
#pull all files

cp -r ../../LayeredGraphAPI .
cp -r ../../PubTator .
cp ../../PhenotypeAPI/PhenotypeCorrelationParser.py ./PubTator/

./buildGraph.sh $1

echo "Cleaning up..."
if [ -d ../../HPO_graph_data ]; then
    rm -rf ../../HPO_graph_data
fi

mv ./HPO_graph_data ../../HPO_graph_data
if [ -z "$1" ]; then
    rm -rf HPO_data_files
fi
rm -rf LayeredGraphAPI
rm -rf PubTator