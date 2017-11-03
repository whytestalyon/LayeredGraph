# LayeredGraph
Methods for building layered graphs and running graph traversal algorithms

## Software Requirements
1. Python - developed and tested using Python 3.6.2.
2. Python packages - numpy, scipy

## LayeredGraph Data Structure
### Purpose
The main purpose of this repository is to store a LayeredGraph (LG) data structure that is stored in LayeredGraph.py.  The LG is a directed graph composed of "layers" of nodes that are connected in different way.
For example, one might have a layer for storing a meal ontology and then a layer for storing an ingredients ontology.  Some meals are related to each other (e.g. a hamburger and a cheeseburger)
and some ingredients will similarly be related (e.g. cucumber and pickle).  Also, there are connections between the layers such as a meal being composed of specific ingredients.  Once
all this information is present in a layered, algorithms can be run over the graph to identify closely related things.  For example, someone might have bread and turkey, the appropriate
graph algorithm might help them realize that adding cheese will create a turkey and cheese sandwich.

### Building a LayeredGraph for Random Walk with Restart (RWR)
1. Create a graph and load all nodes into the graph.
```python
#instantiate a graph
import LayeredGraph
lg = LayeredGraph.LayeredGraph()
#the first value is the layer, second is the node label
lg.addNode('ingredient', 'cheese')
...
lg.addNode('meal', 'sandwich')
#call this only when all nodes have been added
lg.finalizeNodeList()
```
2. Add all edges to the graph.
```python
#values are layer1, node1, layer2, node2, edge weight, and whether the edge is undirected
lg.addEdge('ingredient', 'cheese', 'meal', 'sandwich', 1, True)
...
lg.addEdge('ingredient', 'turkey', 'meal', 'turkey sandwich', 1, True)
lg.addEdge('meal', 'sandwich', 'meal', 'turkey sandwich', 1, True)
```
3. Build the graph transition matrix.  The first line describes how different layers interact during random walk, the second actually generates the matrix necessary for random walk.
```python
#sets transition rates from one layer to another as equal
lg.setGraphJumpEqual()
#calculates the final transition probabilities for a given LayeredGraph
lg.calculateTransitionMatrix()
```
4. Run the RWR algorithm.  This example will return a ranking of all 'meal' nodes in the graph based on the RWR algorithm.
```python
#define the RWR parameters; we have cheese and bread and want to know what meals are closest related to those ingredients
startNodes={('ingredient', 'cheese'):1.0, ('ingredient', 'bread'):1.0}
restartProb=0.1
targetLayers=set(['meal'])
lg.RWR_rank(startNodes, restartProb, targetLayers)
```

## HPO Layered Graph
The LayeredGraphTest.py file contains the information for building an Human Phenotype Ontology (HPO) layered graph.  The intent of the graph is to provide a mechanism for obtaining
genes that are closely related to a particular set of phenotypes.  In this test file, there are several subroutines that are capable of loading the HPO graph into a layer, the genes
from HPO into another layer, and the connections from HPO to genes.  Additionally, there are several other subroutines that load other layers such as disease layers, CPDB, or Omnipath.
In practice, we did not find that any of those layers improved the result, so they are all toggled off in the current code.  This file also includes the methods for running a set of
test cases that were used to compare the results.

Currently, this file is hard-coded to datasets that have been downloaded and processed on a particular machine.  There are several other normalization and weighting decisions built
into the uploaded HPO to gene graph that are not talked about here.  For more information, refer to the source code or email jholt@hudsonalpha.org.

## HPO Server
### Local install
Easy installation can be obtained by running the supplied ```install.sh```.  This will install necessary python3 packages and download
data from the Morgan shared server for use.

Once the pre-requisites are installed, simply run ```python3 application.py``` and navigate to [http://127.0.0.1:5000/search](http://127.0.0.1:5000/search) to access the GUI.  This will return a JSON
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
