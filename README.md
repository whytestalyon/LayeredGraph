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
### Docker install
Assuming you have docker installed on your machine it's fairly simple to build and run the server within docker (no special
setup needed on your machine besides having docker installed).

First you'll need to build the container (make sure docker is running on your machine):
```bash
cd docker/application
./buildApplicationContainer.sh
```

Then you'll need to build the fixture image that contains all of the graph data:
```bash
cd docker/fixtures
./buildFixturesContainer.sh
```

This will first prompt for login info to download data from the Morgan shared server for use. Then it'll call docker to 
build a volume container containing the graph data that will be used by the server.

**NOTE:** If changes to the graph data or the application code need to be reloaded you'll have to rebuild the corresponding 
container first via the respective method listed above and then rerun the application (see below steps on starting the app and stopping it).

### Running the Application
To run the server enter the following command into a terminal (making sure you are in the root directory of the project):
```bash
docker-compose up -d
```
at this point the server will be running! To access it simply get the IP address of your docker machine, if you're on 
non-native docker for Mac (`docker-machine ip dev` or `docker-machine ip default`), or use localhost on other platforms (and native docker on Mac)
and you can access the UI from the URL [http://localhost:5000](http://localhost:5000). 

### Debugging the application
The nice thing about docker is that it performs some nice magic with respect to logging for us. To get a live feed of the 
output from the server on stdout and stderr simply use the following command:
```bash
docker logs -f layeredgraph_rest_1
```
and this will give you a live feed of the output from the server. To stop monitoring the logs simply just hit `ctrl + c`.

### Stopping the server
To stop the server enter the following command into a terminal (making sure you are in the root directory of the project):
```bash
docker-compose down
```
at this point the server will be stopped, and the containers will be cleaned up!