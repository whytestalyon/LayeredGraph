#!/usr/bin/env bash
#temporarilly copy the application code here so it's included in the docker build context
cp -r ../../static .
cp -r ../../templates .
cp -r ../../LayeredGraphAPI .
cp -r ../../application.py .
cp -r ../../PhenotypeAPI .
cp -r ../../PubTator .

#build the application container
docker build -t docker-registry.haib.org/sdi/layered-graph-server .

#remove temporary application code from this directory
rm -rf static
rm -rf templates
rm -rf LayeredGraphAPI
rm -rf PhenotypeAPI
rm -rf PubTator
rm -f application.py
