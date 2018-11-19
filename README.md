# PyxisMap
Methods for building layered graphs and running graph traversal algorithms for analyzing networks in genetics.

## LayeredGraph API
```LayeredGraphAPI/LayeredGraph.py``` contains an API for building a graph containing multiple types or "layers" of nodes (for example, genes and HPO terms).  The
following graphs are all built using this API.  For details on how the API works, refer to the README within the LayeredGraph subdirectory.

### Software Requirements
1. Python - developed and tested using Python 3.6.2.
2. Python packages - numpy, scipy, pronto, networkx
3. Docker and Docker Compose

### Hardware Requirements
1. Docker needs to have at least 4GB of RAM allocated to it from the host machine
2. At minimum Docker needs 2 CPUs allocated to it
3. At minimum DOcker needs 64GB of disk allocated to it

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

### Official Demo Server
We offer an official demo server stood up at [http://pyxis.hudsonalpha.org](http://pyxis.hudsonalpha.org). This is an up-to-date instance of the latest release PyxisMap.
This instance can be used to validate a local deployment of the application.    

### Docker Quick-start
Assuming you have docker installed on your machine it's fairly simple to build and run the server within docker (**IMPORTANT:** Docker needs at least 4GB of memory allocated for building the images).  Here are the instructions for a local build and deployment:

```bash
git clone https://github.com/HudsonAlpha/LayeredGraph.git
cd LayeredGraph/docker/application
./buildApplicationContainer.sh
cd ../fixtures
./buildFixturesContainer.sh
cd ../..
docker-compose up -d
``` 

Navigate to [http://localhost:5000](http://localhost:5000) to see the server running.  More details are available on the [Docker Deployment wiki page](https://github.com/HudsonAlpha/LayeredGraph/wiki/Docker-Deployment).

## License
Please refer to the 'About' section of the PyxisMap UI for information on the licence, or access it directly [here](/static/license/PYXISMAP%20LICENSE%20AND%20TERMS%20OF%20USE.docx)