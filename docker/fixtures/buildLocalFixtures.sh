#!/usr/bin/env bash
#pull all files

cp -r ../../LayeredGraphAPI .
cp -r ../../PubTator .

./buildGraph.sh

echo "Cleaning up..."
if [ -d ../../HPO_graph_data ]; then
    rm -rf ../../HPO_graph_data
fi
mv ./HPO_graph_data ../../HPO_graph_data
rm -rf HPO_data_files
rm -rf LayeredGraphAPI
rm -rf PubTator