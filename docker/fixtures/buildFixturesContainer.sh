#!/usr/bin/env bash

cp -r ../../LayeredGraphAPI .
cp -r ../../PubTator .
cp ../../PhenotypeAPI/PhenotypeCorrelationParser.py ./PubTator/

docker build -t docker-registry.haib.org/sdi/layered-graph-fixtures .
if [ $? -ne 0 ]; then
    echo "Failed to build the fixtures container..."
else
    echo "Fixtures image built!"
fi

rm -rf ./LayeredGraphAPI
rm -rf ./PubTator