# PyxisMap
Methods for building layered graphs and running graph traversal algorithms for analyzing networks in genetics.

## LayeredGraph API
```LayeredGraph/LayeredGraph.py``` contains an API for building a graph containing multiple types or "layers" of nodes (for example, genes and HPO terms).  The
following graphs are all built using this API.  For details on how the API works, refer to the README within the LayeredGraph subdirectory.

### Software Requirements
1. Python - developed and tested using Python 3.6.2.
2. Python packages - numpy, scipy

## Pre-constructed Graphs
### Human Phenotype Ontology (HPO)
The first graph contains two layers: a phenotype layer and a gene layer.  The intent of the graph is to provide a mechanism for obtaining
genes that are closely related to a particular set of phenotypes.  The graph is built using the full HPO database and weights are calculated
between HPO nodes based on the number of shared gene associations.  Phenotype-to-gene edges are built using a combination of HPO annotations and
a scaling weight from PubTator.

### Protein-Protein Interaction (PPI)
The second graph contains a single layer for genes.  The intent of the graph is to provide a mechanism for identifying genes that are interacting with
a set of other genes.  The graph is built using protein-protein interactions provided by ConsensusPathDB where the proteins have been translated into the corresponding
human gene label.  Gene-to-gene edge weights are calculated based on the maximum confidence provided by ConsensusPathDB for the particular relationship.

## HPO Server
### Local install
Easy installation can be obtained by running the supplied ```install.sh```.  This will install necessary python3 packages and download both pre-constructed
graphs from the Morgan shared server for use.

Once the pre-requisites are installed, simply run ```python3 application.py``` and navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000) to access the GUI.  This will return a JSON
output that can also be accessed programmatically via [http://127.0.0.1:5000/rank](http://127.0.0.1:5000/rank).

### Docker install
Assuming you have docker installed on your machine it's fairly simple to build and run the server within docker (no special
setup needed on your machine besides having docker installed).

First you'll need to build the container (make sure docker is running on your machine):
```bash
./buildDockerContainer.sh
```
This will first prompt for login info to download data from the Morgan shared server for use. Then it'll call docker to 
build the container where the server will be running.

To run the server enter the following command into a terminal:
```bash
docker run -d --name lg_server --rm -p 5000:5000 docker-registry.haib.org/sdi/layered-graph-server
```
at this point the server will be running! To access it simply get the IP address of your docker machine, if you're on 
non-native docker for Mac (`docker-machine ip dev` or `docker-machine ip default`), or use localhost (127.0.0.1) on other platforms
and you can access it from the URL listed above (replacing the IP address of 127.0.0.1 with that of the address of your docker machine on Macs).
