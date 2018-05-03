# PyxisMap
Methods for building layered graphs and running graph traversal algorithms for analyzing networks in genetics.

## LayeredGraph API
```LayeredGraph/LayeredGraph.py``` contains an API for building a graph containing multiple types or "layers" of nodes (for example, genes and HPO terms).  The
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

#### Local deployment for live development
To run the server such that you can do local development on the files and have the changes rendered immediately in a terminal
(making sure you are in the root directory of the project):
```bash
docker-compose -f local-docker-compose.yml up
``` 

#### Regular deployment
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


## REST API
To allow for extension and greater use of the application the web application was built on top of a RESTful API.

### Paths
#### Query HPO terms
| Field | Value |
| ----- | ----- |
| Description | returns a JSON presentation of the terms that match user typed text |
| Method | GET |
| Path | `/terms` |
| Return | JSON |
| Query Parameters | `term` = string of text to search for; can be HPO ID or disease/phenotype text |  

Example Call: 

`/terms?term=HP%3A0001234`

Example results: 
```json
{
  "results": [
    {
      "id": "HP:0001234",
      "text": "HP:0001234 Hitchhiker thumb; Hitchhiker thumb; Abducted thumb"
    }
  ]
}
```

#### Query Gene Symbols
| Field | Value |
| ----- | ----- |
| Description | returns a JSON presentation fo the genes that match the supplied search string |
| Method | GET |
| Path | `/genes` |
| Return | JSON |
| Query Parameters | `term` = string of gene symbol to search for |  

Example Call: 

`/genes?term=A2M`

Example results: 
```json
{
  "results": [
    {
      "id": "A2M",
      "text": "A2M"
    },
    {
      "id": "A2M-AS1",
      "text": "A2M-AS1"
    },
    {
      "id": "A2ML1",
      "text": "A2ML1"
    },
    {
      "id": "A2ML1-AS1",
      "text": "A2ML1-AS1"
    },
    {
      "id": "A2ML1-AS2",
      "text": "A2ML1-AS2"
    },
    {
      "id": "A2MP1",
      "text": "A2MP1"
    }
  ]
}
```

#### Query Phenotype to Gene correlation
| Field | Value |
| ----- | ----- |
| Description | Given an input gene and phenotype terms return the publications and phenotypes PyxisMap associates with it |
| Method | POST |
| Path | `/phenotypegene` |
| Return | JSON |
| Payload | JSON |  

Example Payload:

```json
{
  "gene": "SPTBN2",
  "phenotypes": [
    "HP:0001251"
  ]
}
```

Example Call: 

`/phenotypegene`

Example results: 
```json
{
  "data": [
    {
      "count": 1,
      "pmid": "10417284",
      "terms": [
        "HP:0001251"
      ]
    },
    {
      "count": 1,
      "pmid": "10434652",
      "terms": [
        "HP:0001251"
      ]
    },
    {
      "count": 1,
      "pmid": "10522902",
      "terms": [
        "HP:0001251"
      ]
    }
  ]
}
```

#### Query Rank Genes
| Field | Value |
| ----- | ----- |
| Description | Given an input list of HPO IDs produce a ranked list of genes associated with said terms |
| Method | POST |
| Path | `/deeprank` |
| Return | JSON |
| Payload | JSON |  

Example Payload:

```json
["HP:0001251", "HP:0001744"]
```

Example Call: 

`/deeprank`

Example results: 
```json
{
  "rankings": [
    {
      "HP:0001251": [
        444,
        0.0009026369410655975
      ],
      "HP:0001744": [
        1,
        0.0027041699843090034
      ],
      "label": "MVK",
      "nodeType": "gene",
      "rank": 1,
      "weight": 0.0025840677844297085
    },
    {
      "HP:0001251": [
        59,
        0.0012961879743216367
      ],
      "HP:0001744": [
        12,
        0.0023167158518850595
      ],
      "label": "SCYL1",
      "nodeType": "gene",
      "rank": 2,
      "weight": 0.0022486806629565823
    },
    ...
  ],
  "usedTerms": [
    "HP:0001251",
    "HP:0001744"
  ],
  "missingTerms": []
}
```

#### Annotate text with HPO terms
| Field | Value |
| ----- | ----- |
| Description | Given an input of text produce HPO IDs of terms found within the text (interface to NCBO annotator) |
| Method | GET |
| Path | `/text/annotate` |
| Return | JSON |
| Query Parameter | `indications` = the text to be annotated with HPO IDs |  

Example Call: 

`/text/annotate?indications=%22HP%3A0001251%22%3A%20%22Ataxia%3B%20Cerebellar%20ataxia%22%2C%0A%20%20%20%20%22HP%3A0001744%22%3A%20%22Splenomegaly%3B%20Increased%20spleen%20size%22%2C%0A%20%20%20%20%22HP%3A0001945%22%3A%20%22Fever%3B%20Fever%3B%20Hyperthermia%3B%20Pyrexia%22%2C%0A%20%20%20%20%22HP%3A0031796%22%3A%20%22Recurrent%3B%20Intermittent%22`

Example results: 
```json
{
  "annotatorStatus": 200,
  "terms": {
    "HP:0001251": "Ataxia; Cerebellar ataxia",
    "HP:0001744": "Splenomegaly; Increased spleen size",
    "HP:0001945": "Fever; Fever; Hyperthermia; Pyrexia",
    "HP:0031796": "Recurrent; Intermittent"
  }
}
```

## License

This software is free for use for academic and non-profit endeavors under the GPL v3. Please contact bwilk@hudsonalpha.org or jholt@hudsonalpha.org to inquire about commercial use. 
