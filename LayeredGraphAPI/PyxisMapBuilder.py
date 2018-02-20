#!/usr/bin/env python3
'''
Usage: python3 PyxisMapBuilder.py

This file contains the methods to construct the LayeredGraph used by PyxisMap.  Pathing is currently local to Matt Holt's desktop.
'''

import json
import numpy as np
import pickle

import LayeredGraph

def createPyxisMapGraph(hpoPhenoToGenoFN, graphStructureFN):
    '''
    This function will create an the graph we use for HPO term to gene
    @param hpoPhenoToGenoFN - the file containing relationships between phenotypes and genotype; typically ./HPO_dl/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
    @param graphStructureFN - the OBO file from HPO; typically ./HPO_dl/hp.obo
    @return - an instance of LayeredGraph containing the constructed graph
    '''
    #first, load in the phenotype to gene mappings
    print('Loading HPO to gene information...')
    p2g = loadPhenoToGeno(hpoPhenoToGenoFN)
    print(len(p2g), 'p2g associations')
    
    #second, build a graph
    nodes, edges, altIDMap = loadGraphStructure(graphStructureFN)
    
    #make sure anything in p2g is stored as the main ID
    for k in p2g:
        assert(k not in altIDMap)
    
    #make sure all our edges only use main IDs also
    for source in edges:
        assert(source not in altIDMap)
        for dest in edges[source]:
            assert(dest not in altIDMap)
    
    #push all p2g info up the HPO graph
    p2g = pushP2gUp(p2g, nodes, edges)
    
    #pubtator relations are stored as (entrez ID, HPO term) pairs, so get the entrez ID dictionary and then get the weights for use downstream
    ed, aliases = getEntrezDict()
    p2gWeights = getP2GWeights()
    
    #init
    mg = LayeredGraph.LayeredGraph()
    
    #add all of our nodes to the graph
    print('Adding HPO nodes...')
    for n in nodes:
        mg.addNode('HPO', n)
    
    print('Gathering gene nodes from HPO...')
    geneSet = set([])
    for p in p2g:
        geneSet |= p2g[p]
    
    print('Gathering gene nodes from Pubtator...')
    pubGeneSet = set([])
    for eid, hpo in p2gWeights:
        if (eid in ed):
            pubGeneSet.add(ed[eid])
    
    #combine the HPO gene set and the pubtator gene set then add those
    print('Adding gene nodes to the graph...')
    comboSet = geneSet | pubGeneSet
    print(str(len(comboSet))+' genes identified for the graph...')
    for g in comboSet:
        mg.addNode('gene', g)
    
    #all nodes are in now
    mg.finalizeNodeList()
    
    #add p2p edges; bidirectional and weighted based on gene overlaps from HPO
    print('Adding HPO edges...')
    for parent in list(edges.keys()):
        for child in edges[parent]:
            #scale with gene similarity
            similarity = len(p2g.get(parent, set([])) & p2g.get(child, set([])))
            #similarity = np.log2(len(p2g.get(parent, set([])) & p2g.get(child, set([]))) + 1)
            mg.addEdge('HPO', parent, 'HPO', child, similarity+1, True)
    
    #all p2g edges are directional
    print('Adding p2g edges...')
    print('WARNING: combined ignoring flags')
    for p in p2g:
        for g in p2g[p]:
            #first, add a single constant weight to each HPO to gene relationship of 1.0
            mg.addEdge('HPO', p, 'gene', g, 1, False)
    
    #identify the Pubtator pairing with the largest weight
    maxValue = 0.0
    for (eid, hpo) in p2gWeights:
        if (hpo not in nodes):
            alias = aliases.get(hpo, None)
            if alias == None:
                print('Warning: '+hpo+' not found in HPO database.')
                continue
            else:
                print('Swapping '+hpo+' for alias '+alias)
                hpo = alias
        
        if p2gWeights[(eid, hpo)] > maxValue:
            maxValue = p2gWeights[(eid, hpo)]
    
    #update any edges that have Pubtator annotations with a weighted fraction from [0.0, 1.0]
    for (eid, hpo) in p2gWeights:
        if (hpo not in nodes):
            alias = aliases.get(hpo, None)
            if alias == None:
                print('Warning: '+hpo+' not found in HPO database.')
                continue
            else:
                print('Swapping '+hpo+' for alias '+alias)
                hpo = alias
        
        #if this has an HPO to gene relationship, add 1.0 to the weight (since its a replacement for the previous weight)
        if (eid in ed):
            g = ed[eid]
            if (g in p2g.get(hpo, [])):
                mg.addEdge('HPO', hpo, 'gene', g, 1.0+(p2gWeights[(eid, hpo)] / maxValue), False)
            else:
                mg.addEdge('HPO', hpo, 'gene', g, p2gWeights[(eid, hpo)] / maxValue, False)
        
    #now finish out everything
    print('Setting graph jump equal...')
    mg.setGraphJumpEqual()
    
    print('Calculating final transition matrix...')
    mg.calculateTransitionMatrix()
    
    return mg

def loadPhenoToGeno(fn):
    '''
    This function parses the HPO phenotype file and returns HPO->gene information
    @param fn - the filename for the HPO file
    @return - a dictionary where the key is an HPO term and the value is a set of gene names
    '''
    #return a dictionary of p2g
    ret = {}
    fp = open(fn, 'r')
    
    #capture the header then read terms
    fp.readline()
    for l in fp:
        pieces = l.strip('\n').split('\t')
        hpo = pieces[0]
        geneName = pieces[3].upper()
        
        if hpo not in ret:
            ret[hpo] = set([geneName])
        else:
            ret[hpo].add(geneName)
    fp.close()
    
    return ret

def loadGraphStructure(fn):
    '''
    This function parses the HPO graph structure and returns information to build the corresponding graph
    @param fn - the .obo file containing the HPO structure
    @return - tuple (nodes, edges, altIdToMain)
        nodes - the set of HPO terms that are nodes
        edges - a dictionary where key is a parent HPO term and value is a set of child HPO terms
        altIdToMain - a dictionary where key is an old HPO term and value is the current HPO term (for cases where terms were consolidated)
    '''
    fp = open(fn, 'r')
    
    nodes = set([])
    edges = {}
    altIdToMain = {}
    
    #TODO: unhandled: alternate ID's
    inTerm = False
    for l in fp:
        if inTerm:
            if l == '\n':
                nodeID = None
                inTerm = False
            elif l[0:4] == 'id: ':
                nodeID = l[4:].strip('\n')
                nodes.add(nodeID)
            elif l[0:6] == 'is_a: ':
                assert(nodeID != None)
                parentID = l[6:16]
                if parentID in edges:
                    edges[parentID].add(nodeID)
                else:
                    edges[parentID] = set([nodeID])
            elif l[0:8] == 'alt_id: ':
                assert(nodeID != None)
                altID = l[8:18]
                altIdToMain[altID] = nodeID
            else:
                pass
            
        elif l[0:6] == '[Term]':
            inTerm = True
    
    fp.close()
    
    return nodes, edges, altIdToMain

def pushP2gUp(p2g, nodes, edges):
    '''
    This function push the HPO phenotype-to-gene information completely up the "tree"-like graph structure
    @param p2g - the dictionary of HPO terms to genes
    @param nodes - the set of nodes in the HPO graph
    @param edges - the set of edges where key is the parent node in the HPO graph and the value is a set of children
    @return - a modified p2g that pushes all phenotype info full up the tree
    '''
    print('Running pushP2gUp(...)')
    analyzed = set([])
    rootNode = 'HP:0000001'
    stack = [rootNode]
    
    while len(stack) > 0:
        currNode = stack[-1]
        for end in edges.get(currNode, set([])):
            if end not in analyzed:
                stack.append(end)
        
        if currNode == stack[-1]:
            #nothing was added
            x = len(p2g.get(currNode, set([])))
            for end in edges.get(currNode, set([])):
                p2g[currNode] = p2g.get(currNode, set([])) | p2g.get(end, set([]))
            stack.pop()
            analyzed.add(currNode)
    
    return p2g

def getEntrezDict():
    '''
    This function created dictionaries that go from Entrez ID to gene name and vice-versa
    @return tuple (ed, aliases)
        ed - key is entrez ID, value is the main geneName
        aliases - key is an entrez ID, value is a list of aliases
    '''
    ed = {}
    aliases = {}
    fp = open('./HPO_data_files/hgnc_complete_set.txt', 'r')
    
    headers = fp.readline().rstrip()
    hDict = {l : i for i, l in enumerate(headers.split('\t'))}
    geneIndex = hDict['symbol']
    alIndex = hDict['alias_symbol']
    entrezIndex = hDict['entrez_id']
    
    for l in fp:
        pieces = l.rstrip().split('\t')
        geneName = pieces[geneIndex].strip('"').upper()
        if len(pieces) < 19:
            continue
        eid = pieces[entrezIndex]
        ed[eid] = geneName
        al = pieces[alIndex].strip('"')
        if al != '':
            aliases[eid] = al.split('|')
    fp.close()
    
    return ed, aliases

def getP2GWeights():
    '''
    This will calculate the p2g weights based on pubtator stuff
    @return - dictionary where key is (entrez ID, hpo term) and value is the number of pmids attached to that pair
    '''
    p2gSets = {}
    fp = open('./HPO_data_files/gene2phenotype.json', 'r')
    
    for l in fp:
        j = json.loads(l)
        gene = j['geneId'].upper()
        hpoTerm = j['hpoId']
        pmids = j['pmids']
        p2gSets[(gene, hpoTerm)] = len(pmids)
    
    fp.close()
    return p2gSets

def createHPOWeights(hpoPhenoToGenoFN, graphStructureFN):
    '''
    This function will create weights for each HPO term
    @param hpoPhenoToGenoFN - the file containing relationships between phenotypes and genotype; typically ./HPO_dl/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
    @param graphStructureFN - the OBO file from HPO; typically ./HPO_dl/hp.obo
    @return - a dictionary where key is an HPO term and value is the weight for that term based on the calculation
    '''
    #first, load in the phenotype to gene mappings
    print('Loading HPO to gene information...')
    p2g = loadPhenoToGeno(hpoPhenoToGenoFN)
    print(len(p2g), 'p2g associations')
    
    #second, build a graph
    nodes, edges, altIDMap = loadGraphStructure(graphStructureFN)
    
    #make sure anything in p2g is stored as the main ID
    for k in p2g:
        assert(k not in altIDMap)
    
    #make sure all our edges only use main IDs also
    for source in edges:
        assert(source not in altIDMap)
        for dest in edges[source]:
            assert(dest not in altIDMap)
    
    #push all p2g info up
    p2g = pushP2gUp(p2g, nodes, edges)
    
    #HPO graph is fully loaded, now we can calculate weights
    hpoScores = calculateHpoScores(nodes, edges)
    return hpoScores, p2g

def calculateHpoScores(nodes, edges):
    '''
    This function will actually calculate the HPO weights based on their placement in the HPO ontology.  In general, terms closer to the root will have a lower weight
    than those near the leaves.  This is intended to reflect how some terms are more specific than others and should be given higher weights.
    @param nodes - the nodes (aka, HPO terms) from the HPO ontology
    @param edges - the edges connecting HPO terms in the ontology
    @return - a dictionary where key is an HPO term and value is the weight for that term based on the calculation
    '''
    #storage for the number of leaves and subsumers
    leafDict = {}
    subsumerDict = {}
    
    #initialize the stack to the term 'all'
    analyzed = set([])
    rootNode = 'HP:0000001'
    stack = [rootNode]
    
    #while the stack isn't empty
    while len(stack) > 0:
        #look at the last node and check if all of its kids are finished computing
        currNode = stack[-1]
        for end in edges.get(currNode, set([])):
            if end not in analyzed:
                stack.append(end)
        
        if currNode == stack[-1]:
            #nothing was added
            leaves = set([])
            if len(edges.get(currNode, set([]))) == 0:
                #no children, so it is a leaf node
                leaves.add(currNode)
            else:
                #it has children, collect all of its children's leaves
                for end in edges.get(currNode, set([])):
                    leaves = leaves | leafDict[end]
            
            #store, pop, and mark as analyzed
            leafDict[currNode] = leaves
            stack.pop()
            analyzed.add(currNode)
    
    #calculate the parents of each node
    parents = {}
    for source in edges:
        for dest in edges[source]:
            if dest not in parents:
                parents[dest] = set([])
            parents[dest].add(source)
    
    #now go through and calculate the number of subsumers of each node
    analyzed = set([])
    queue = sorted(nodes)
    while len(queue) > 0:
        currNode = queue[0]
        for par in parents.get(currNode, set([])):
            if par not in analyzed:
                queue.append(queue.pop(0))
                break
        
        if currNode == queue[0]:
            #all of this node's parents have been analyzed, we can calculate it
            subsume = set([currNode])
            for parent in parents.get(currNode, set([])):
                subsume = subsume | parents.get(parent, set([]))
            
            #store, pop, and mark as analyzed
            subsumerDict[currNode] = subsume
            queue.pop(0)
            analyzed.add(currNode)
    
    #finally, do the actual math
    leafCounts = [len(leafDict.get(k, [k])) for k in sorted(nodes)]
    subsumCounts = [len(subsumerDict[k]) for k in sorted(nodes)]
    maxLeaves = max(leafCounts)
    
    #originally I was using these weights which applied e^x 
    hpoScores = {k: np.exp(np.log(subsumCounts[x])+np.log(maxLeaves)-np.log(leafCounts[x])) for x, k in enumerate(sorted(nodes))}
    #hpoScores = {k: (np.log(subsumCounts[x])+np.log(maxLeaves)-np.log(leafCounts[x])) for x, k in enumerate(sorted(nodes))}
    return hpoScores

if __name__ == '__main__':
    #static files we will be using
    hpoPhenoToGenoFN = './HPO_data_files/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'
    graphStructureFN = './HPO_data_files/hp.obo'
    
    #output files
    pickleGraphFN = './HPO_graph_data/multigraph.pickle'
    pickleHPOWeightFN = './HPO_graph_data/multiHpoWeight_biogrid_pushup.pickle'
    pickleP2gFN = './HPO_graph_data/multiP2G_biogrid_pushup.pickle'
    
    #create the graph
    mg = createPyxisMapGraph(hpoPhenoToGenoFN, graphStructureFN)
    
    #load or generate the graph
    print('Pickling PyxisMap...')
    fp = open(pickleGraphFN, 'wb+')
    pickle.dump(mg, fp)
    fp.close()
    
    #load or generate the HPO weights
    print('Generating HPO weights...')
    hpoWeights, p2g = createHPOWeights(hpoPhenoToGenoFN, graphStructureFN)
    
    print('Pickling HPO weights...')
    fp = open(pickleHPOWeightFN, 'wb+')
    pickle.dump(hpoWeights, fp)
    fp.close()
    
    fp = open(pickleP2gFN, 'wb+')
    pickle.dump(p2g, fp)
    fp.close()
    