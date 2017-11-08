
import csv
import glob
import json
import numpy as np
import os
#import pickle
import pickle as pickle
import time

import LayeredGraph

def createNewLayeredGraph(hpoPhenoToGenoFN, graphStructureFN, includeDisease, includeCPDB, includeOmnipath):
    '''
    This function will create an the graph we use for HPO term to gene
    @param hpoPhenoToGenoFN - the file containing relationships between phenotypes and genotype; typically ./HPO_dl/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
    @param graphStructureFN - the OBO file from HPO; typically ./HPO_dl/hp.obo
    @param includeDisease - if True, add disease layer to graph and connect it to phenotypes (from HPO)
    @param includeCPDB - if True, add pathway layer to graph and connect it to genes (from CPDB)
    @param includeOmnipath - if True, add gene-to-gene interactions based on protein-to-gene relationships (from Omnipath)
    @return - an instance of LayeredGraph containing the constructed graph
    '''
    #parameters for how the graph generation should be handled
    PUSHUP = True
    WEIGHTED_P2G = True
    WEIGHTED_P2P = True
    
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
    if PUSHUP:
        p2g = pushP2gUp(p2g, nodes, edges)
    
    if includeDisease:
        d2p = loadDiseaseToPheno()
        print(len(d2p), 'd2p interactions loaded')
    else:
        d2p = []
    
    if includeCPDB:
        cpdb = loadCPDB()
        print(len(cpdb), 'pathways loaded')
    else:
        cpdb = []
        
    if includeOmnipath:
        prots, protEdges, prot2gene = loadOmnipath()
        print(len(protEdges), 'omnipath interactions loaded')
    else:
        prots, protEdges, prot2gene = set([]), set([]), {}
    
    if WEIGHTED_P2G:
        ed, aliases = getEntrezDict()
        p2gWeights = getP2GWeights()
    
    #init
    mg = LayeredGraph.LayeredGraph()
    
    #add all of our nodes to the graph
    print('Adding HPO nodes...')
    for n in nodes:
        mg.addNode('HPO', n)
    
    print('Adding gene nodes...')
    geneSet = set([])
    for p in p2g:
        geneSet |= p2g[p]
    
    if WEIGHTED_P2G:
        print('Using Pubtator-based p2g...')
        pubGeneSet = set([])
        for eid, hpo in p2gWeights:
            pubGeneSet.add(ed[eid])
    
    if True:
        #############
        comboSet = geneSet | pubGeneSet
        print(str(len(comboSet))+' genes identified for the graph...')
        for g in comboSet:
            mg.addNode('gene', g)
        #############
    
    else:
        print(str(len(geneSet))+' genes identified for the graph...')
        for g in geneSet:
            mg.addNode('gene', g)
    
    if includeDisease:
        print('Adding disease nodes...')
        #here is where we would add disease nodes
        for d, p in d2p:
            mg.addNode('disease', d)
    
    if includeCPDB:
        print('Adding pathway nodes...')
        for pathwayID, geneList in cpdb:
            mg.addNode('pathway', pathwayID)
    
    if includeOmnipath:
        print('Adding protein nodes...')
        for prot in prots:
            mg.addNode('protein', prot)
        
        for prot in list(prot2gene.keys()):
            mg.addNode('protein', prot)
    
    #all nodes are in now
    mg.finalizeNodeList()
    
    #add all of our edges to the graph
    print('Adding HPO edges...')
    for parent in list(edges.keys()):
        for child in edges[parent]:
            #first method, all edges have constant weight
            if not WEIGHTED_P2P:
                mg.addEdge('HPO', parent, 'HPO', child, 1, True)
            else:
                #second method, scales with gene similarity
                similarity = len(p2g.get(parent, set([])) & p2g.get(child, set([])))
                #similarity = np.log2(len(p2g.get(parent, set([])) & p2g.get(child, set([]))) + 1)
                mg.addEdge('HPO', parent, 'HPO', child, similarity+1, True)
    
    print('Adding p2g edges...')
    '''
    if WEIGHTED_P2G:
        for (eid, hpo) in p2gWeights:
            if (hpo not in nodes):
                alias = aliases.get(hpo, None)
                if alias == None:
                    print('Warning: '+hpo+' not found in HPO database.')
                    continue
                else:
                    print('Swapping '+hpo+' for alias '+alias)
                    hpo = alias
            
            g = ed[eid]
            if (g in geneSet):
                mg.addEdge('HPO', hpo, 'gene', g, p2gWeights[(eid, hpo)], False)
            else:
                print('Ignoring '+hpo+'-'+g)
    else:
        for p in p2g:
            for g in p2g[p]:
                #first method, all edges have constant weight
                #mg.addEdge('HPO', p, 'gene', g, 1, False)
                
                #second method, scales with gene similarity
                mg.addEdge('HPO', p, 'gene', g, 2, False)
    '''
    
    ################################################################
    print('WARNING: combined ignoring flags')
    for p in p2g:
        for g in p2g[p]:
            #first method, all edges have constant weight
            mg.addEdge('HPO', p, 'gene', g, 1, False)
            
            #second method, scales with gene similarity
            #mg.addEdge('HPO', p, 'gene', g, 2, False)
    
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
        #if np.log2(p2gWeights[(eid, hpo)]) > maxValue:
            maxValue = p2gWeights[(eid, hpo)]
            #maxValue = np.log2(p2gWeights[(eid, hpo)])
        
    for (eid, hpo) in p2gWeights:
        if (hpo not in nodes):
            alias = aliases.get(hpo, None)
            if alias == None:
                print('Warning: '+hpo+' not found in HPO database.')
                continue
            else:
                print('Swapping '+hpo+' for alias '+alias)
                hpo = alias
        
        g = ed[eid]
        if True or (g in geneSet):
            if (g in p2g.get(hpo, [])):
                mg.addEdge('HPO', hpo, 'gene', g, 1.0+(p2gWeights[(eid, hpo)] / maxValue), False)
                #mg.addEdge('HPO', hpo, 'gene', g, 1.0+(np.log2(p2gWeights[(eid, hpo)]) / maxValue), False)
            else:
                mg.addEdge('HPO', hpo, 'gene', g, p2gWeights[(eid, hpo)] / maxValue, False)
                #mg.addEdge('HPO', hpo, 'gene', g, np.log2(p2gWeights[(eid, hpo)]) / maxValue, False)
        else:
            print('Ignoring '+hpo+'-'+g)
    ################################################################
    
    if includeDisease:
        print('Adding d2p edges...')
        for d, p in d2p:
            mg.addEdge('disease', d, 'HPO', p, 1, True)
    
    if includeCPDB:
        print('Adding self edges on genes...')
        for g in geneSet:
            mg.addEdge('gene', g, 'gene', g, 1, False)
            
        print('Adding cpdb edges...')
        for pathwayID, geneList in cpdb:
            for g in geneList:
                if g in geneSet:
                    mg.addEdge('gene', g, 'pathway', pathwayID, 1, True)
    
    if includeOmnipath:
        print('Adding self edges on genes...')
        for g in geneSet:
            mg.addEdge('gene', g, 'gene', g, 1, False)
        
        print('Adding gene to protein edges...')
        for prot in prot2gene:
            for g in prot2gene[prot]:
                if g in geneSet:
                    mg.addEdge('gene', g, 'protein', prot, 1, True)
        
        print('Adding PPI...')
        for source, dest in protEdges:
            mg.addEdge('protein', source, 'protein', dest, 1, False)
    
    #now finish out everything
    #print 'Setting graph jump equal...'
    mg.setGraphJumpEqual()
    
    print('Calculating final transition matrix...')
    mg.calculateTransitionMatrix()
    
    return mg

def createHPOWeights(hpoPhenoToGenoFN, graphStructureFN):
    '''
    This function will create weights for each HPO term
    @param hpoPhenoToGenoFN - the file containing relationships between phenotypes and genotype; typically ./HPO_dl/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
    @param graphStructureFN - the OBO file from HPO; typically ./HPO_dl/hp.obo
    @return - a dictionary where key is an HPO term and value is the weight for that term based on the calculation
    '''
    #parameters for how the graph generation should be handled
    PUSHUP = True
    
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
    if PUSHUP:
        p2g = pushP2gUp(p2g, nodes, edges)
    
    #HPO graph is fully loaded, now we can calculate weights
    hpoScores = calculateHpoScores(nodes, edges)
    return hpoScores, p2g

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
        geneName = pieces[3]
        
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
                #if len(nodes) >= 10:
                #    break
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
                #print l.strip('\n')
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

def deriveParents(nodes, edges):
    '''
    This function will calculate all parent terms for each HPO node
    @param nodes - the set of nodes in the HPO graph
    @param edges - the set of edges where key in the parent node in the HPO graph and the value is a set of children
    @return - a dictionary where key is an HPO term and value is a set of nodes that are ancestors of that term
    '''
    #first calculate all immediate parents
    allKeys = set([])
    revEdges = {}
    for parent in edges:
        allKeys.add(parent)
        for child in edges[parent]:
            if (child not in revEdges):
                revEdges[child] = set([])
            revEdges[child].add(parent)
            allKeys.add(child)
    
    analyzed = set([])
    topush = list(allKeys)
    ret = {}
    
    while len(topush) > 0:
        currNode = topush.pop(0)
        if revEdges.get(currNode, set([])).issubset(analyzed):
            #all parent terms have finished
            analyzed.add(currNode)
            ret[currNode] = set([])
            for parent in revEdges.get(currNode, set([])):
                ret[currNode] |= ret[parent]
        else:
            #still waiting on parent terms
            topush.append(currNode)
    
    return ret

def calculateHpoScores(nodes, edges):
    '''
    This function will actually calculate the HPO weights based on their placement in the HPO ontology.  In general, terms closer to the root will have a lower weight
    than those near the leaves.  This is intended to reflect how some terms are more specific than others and should be given higher weights.
    @param nodes - the nodes (aka, HPO terms) from the HPO ontology
    @param edges - the edges connecting HPO terms in the ontology
    @return - a dictionary where key is an HPO term and value is the weight for that term based on the calculation
    '''
    leafDict = {}
    subsumerDict = {}
    
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
            #x = len(p2g.get(currNode, set([])))
            #for end in edges.get(currNode, set([])):
            #    p2g[currNode] = p2g.get(currNode, set([])) | p2g.get(end, set([]))
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
    
    parents = {}
    for source in edges:
        for dest in edges[source]:
            if dest not in parents:
                parents[dest] = set([])
            parents[dest].add(source)
    
    analyzed = set([])
    queue = sorted(nodes)
    while len(queue) > 0:
        currNode = queue[0]
        for par in parents.get(currNode, set([])):
            if par not in analyzed:
                queue.append(queue.pop(0))
                break
        
        if currNode == queue[0]:
            subsume = set([currNode])
            for parent in parents.get(currNode, set([])):
                subsume = subsume | parents.get(parent, set([]))
            
            #store, pop, and mark as analyzed
            subsumerDict[currNode] = subsume
            queue.pop(0)
            analyzed.add(currNode)
    
    leafCounts = [len(leafDict.get(k, [k])) for k in sorted(nodes)]
    subsumCounts = [len(subsumerDict[k]) for k in sorted(nodes)]
    maxLeaves = max(leafCounts)
    
    #originally I was using these weights which applied e^x 
    hpoScores = {k: np.exp(np.log(subsumCounts[x])+np.log(maxLeaves)-np.log(leafCounts[x])) for x, k in enumerate(sorted(nodes))}
    #hpoScores = {k: (np.log(subsumCounts[x])+np.log(maxLeaves)-np.log(leafCounts[x])) for x, k in enumerate(sorted(nodes))}
    return hpoScores

def loadDiseaseToPheno():
    '''
    This function reads the HPO disease-phenotype file and returns a list of disease-phenotype pairs
    @return - list of tuples (disease, HPO)
    '''
    ret = []
    fp = open('/Users/matt/data/HPO_dl/phenotype_annotation.tab', 'r')
    
    allowedDBs = set(['OMIM'])
    #allowedDBs = set(['DECIPHER'])
    #allowedDBs = set(['ORPHA'])
    #allowedDBs = set(['OMIM', 'DECIPER', 'ORPHA'])
    
    for l in fp:
        pieces = l.rstrip().split('\t')
        if pieces[0] in allowedDBs:
            disease = pieces[5]
            pheno = pieces[4]
            
            ret.append((disease, pheno))
            
    fp.close()
    return ret

def loadCPDB():
    '''
    This function reads the CPDB file and returns a list of pathway-gene pairs
    @return - list of tuples (pathway label, set of genes)
    '''
    ret = []
    fp = open('/Users/matt/Downloads/CPDB_pathways_genes-3.tab', 'r')
    fp.readline()
    
    allowedSources = set(['KEGG'])
    sources = set([])
    for l in fp:
        pieces = l.rstrip().split('\t')
        pathway, external_id, source, geneList = pieces
        geneList = set(geneList.split(','))
        if source in allowedSources:
            ret.append((external_id+'_'+source, geneList))
            sources.add(source)
        
    fp.close()
    
    print('sources: '+str(sources))
    return ret

def loadOmnipath():
    '''
    This loads the Omnipath datasets which are protein protein interactions (or could be used as g2g interactions also)
    @return tuple (prots, protEdges, prot2gene)
        prots - the protein node lists
        protEdges - the protein to protein edge list
        prot2gene - the connections between proteins and genes
    '''
    #this method did not seem to actually help, maybe try an alternate g2g interaction network?
    #print 'Loading gene to gene information from omnipath...'
    proteinToEnsFN = '/Users/matt/data/HUMAN_9606_idmapping_selected.tab'
    omnipathFN = '/Users/matt/data/interactions_omnipath.tsv'
    prots, protEdges = loadProteinInteractions(omnipathFN)
    #print len(prots), len(protEdges)
    
    p2e = loadProteinToEns(proteinToEnsFN)
    #print len(p2e)
    
    ens2GeneFN = '/Users/matt/data/ens2gene.txt'
    e2g = loadEns2Gene(ens2GeneFN)
    #print len(e2g)
    
    prot2gene = condenseProt2Gene(p2e, e2g)
    #print len(prot2gene)
    
    return prots, protEdges, prot2gene

def loadProteinInteractions(fn):
    '''
    This parses the actual omnipath interaction file
    fn - the path to the file
    @return tuple (prots, edges)
        prots - a list of proteins
        edges - a list of edges
    '''
    fp = open(fn, 'r')
    
    prots = set([])
    edges = set([])
    
    fp.readline()
    for l in fp:
        pieces = l.rstrip().split('\t')
        source = pieces[0]
        target = pieces[1]
        is_directed = int(pieces[2])
        
        prots.add(source)
        prots.add(target)
        edges.add((source, target))
        if is_directed == 0:
            edges.add((target, source))
    
    fp.close()
    
    return (prots, edges)

def loadEns2Gene(fn):
    '''
    This parses the file connecting ensemblIDs to gene names
    @param fn - the file path
    @return - a dictionary from ensemblID to gene name
    '''
    ret = {}
    fp = open(fn, 'r')
    
    fp.readline()
    for l in fp:
        pieces = l.rstrip().split('\t')
        ens = pieces[0]
        geneName = pieces[1]
        ret[ens] = geneName
    
    fp.close()
    
    return ret

def loadProteinToEns(fn):
    ret = {}
    fp = open(fn, 'r')
    
    for l in fp:
        pieces = l.rstrip().split('\t')
        #print pieces
        if len(pieces) > 18 and pieces[18] != '':
            proteinLabel = pieces[0]
            ensLabels = pieces[18].split('; ')
            
            for la in ensLabels:
                #print proteinLabel, la
                assert(len(la) == 15)
            assert(proteinLabel not in ret)
            
            ret[proteinLabel] = ensLabels
            
    fp.close()
    
    return ret

def getP2GWeights():
    '''
    This will calculate the p2g weights based on pubtator stuff
    @return - dictionary where key is (entrez ID, hpo term) and value is the number of pmids attached to that pair
    '''
    p2gSets = {}
    fp = open('/Users/matt/data/HPO_dl/gene2phenotype.json', 'r')
    
    for l in fp:
        j = json.loads(l)
        gene = j['geneId']
        hpoTerm = j['hpoId']
        pmids = j['pmids']
        p2gSets[(gene, hpoTerm)] = len(pmids)
    
    fp.close()
    return p2gSets

def getEntrezDict():
    '''
    This function created dictionaries that go from Entrez ID to gene name and vice-versa
    @return tuple (ed, aliases)
        ed - key is entrez ID, value is the main geneName
        aliases - key is an entrez ID, value is a list of aliases
    '''
    ed = {}
    aliases = {}
    fp = open('/Users/matt/data/HPO_dl/hgnc.txt', 'r')
    fp.readline()
    for l in fp:
        pieces = l.rstrip().split('\t')
        geneName = pieces[1].strip('"')
        if len(pieces) < 19:
            continue
        eid = pieces[18]
        ed[eid] = geneName
        al = pieces[8].strip('"')
        if al != '':
            aliases[eid] = al.split('|')
        #ed[geneName] = eid
        #aliases = pieces[8].strip('"').split('|')
        #for al in aliases:
        #    if al != '':
        #        ed[al] = eid
    fp.close()
    return ed, aliases

def condenseProt2Gene(p2e, e2g):
    ret = {}
    for prot in p2e:
        ensList = p2e[prot]
        
        geneSet = set([])
        for ens in ensList:
            if ens in e2g:
                geneSet.add(e2g[ens])
        
        if len(geneSet) > 0:
            ret[prot] = list(geneSet)
    return ret

def runTests(mg):
    #AFF4:
    #hpoTerms = ['HP:0003196', 'HP:0011951', 'HP:0000280', 'HP:0000518', 'HP:0000527', 'HP:0000648']
    
    cases = []
    
    #Camille Case #1 - ADA, CD3E, RAG1
    hpoTerms = ['HP:0005352', 'HP:0011343', 'HP:0001999', 'HP:0010609', 'HP:0002779']
    geneSearch = ['ADA', 'CD3E', 'RAG1']
    cases.append((hpoTerms, geneSearch))
    
    #Camille Case #2 - PARK2, ATP2B3
    hpoTerms = ['HP:0007272', 'HP:0001300', 'HP:0001317', 'HP:0007240']
    geneSearch = ['PARK2']
    cases.append((hpoTerms, geneSearch))
    
    #Camille Case #3 - AMBN, BFSP2
    hpoTerms = ['HP:0001118', 'HP:0011474', 'HP:0001251', 'HP:0006297', 'HP:0001249', 'HP:0001999', 'HP:0001388']
    geneSearch = ['AMBN', 'BFSP2']
    cases.append((hpoTerms, geneSearch))
    
    #Camille Case #4 - MPO, PGAP2
    hpoTerms = ['HP:0001987', 'HP:0002579']
    geneSearch = ['MPO', 'PGAP2']
    cases.append((hpoTerms, geneSearch))
    
    #load in these cases as well
    cases2 = loadCasesFromMana()
    cases += cases2
    
    #cases = loadCasesFromMel()
    #exit()
    
    hpoWeights = pickle.load(open('/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle', 'r'))
    
    ranks = []
    for hpoTerms, geneSearch in cases:
        if len(hpoTerms) == 0 or len(geneSearch) == 0:
            continue
        
        print(hpoTerms)
        #startProbs = {('HPO', h) : 1.0 for h in hpoTerms}
        startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
        #print startProbs
        restartProb = 0.1
        rankTypes = set(['gene'])
        rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes)
        
        #second, check the rank of genes in the diffusion
        for gene in geneSearch:
            for x in range(0, len(rankedGenes)):
                if rankedGenes[x][2] == gene:
                    w = rankedGenes[x][0]
                    r = x
                    break
            
            #expand in the event of equals
            l = r
            h = r
            while l > 0 and rankedGenes[l-1][0] == w:
                l -= 1
            while h < len(rankedGenes)-1 and rankedGenes[h+1][0] == w:
                h += 1
            finalRank = (l+h)/2
                
            ranks.append(finalRank)
            print(gene+' ranked '+str(finalRank))
        print()
        
    print('Ranks\t'+str(ranks))
    print('Total_rank\t'+str(np.sum(ranks)))
    print('Mean_rank\t'+str(np.mean(ranks))+' ('+str(np.mean(np.log(np.array(ranks)+1)))+')')
    print('Min_rank\t'+str(np.min(ranks)))
    print('Max_rank\t'+str(np.max(ranks)))
    print()

def loadCasesFromMana():
    HPOFiles = glob.glob('/Users/matt/Downloads/4matt-udn-cases/*/HPO_input.txt')
    reportedGenes = '/Users/matt/Downloads/4matt-udn-cases/UDN_reported_results.csv'
    
    caseToHPO = {}
    
    for fn in HPOFiles:
        caseLabel = fn.split('/')[-2]
        fp = open(fn, 'r')
        hpo = set([])
        for l in fp:
            if l[0:3] == 'HP:':
                hpo.add(l.strip('\n'))
        fp.close()
        
        caseToHPO[caseLabel] = hpo
    
    caseToGenes = {}
    fp = open(reportedGenes, 'r')
    fp.readline()
    for l in fp:
        pieces = l.strip('\n').strip('\r').split(',')
        geneLabel = pieces[0]
        caseLabel = pieces[1]
        classification = pieces[11]
        
        if classification == 'primary':
            if caseLabel in caseToGenes:
                caseToGenes[caseLabel].add(geneLabel)
            else:
                caseToGenes[caseLabel] = set([geneLabel])
    fp.close()
    
    retCases = []
    for k in list(caseToHPO.keys()):
        if k in caseToGenes:
            retCases.append((list(caseToHPO[k]), list(caseToGenes[k])))
    
    return retCases

def testRankings(mg):
    cases = rankCasesFromMana()
    
    hpoWeights = pickle.load(open('/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle', 'rb'))
    #hpoFreq = getHPODiseaseFreq()
    #maxFreq = max([hpoFreq[k] for k in hpoFreq])
    #hpoWeights = {k : 1.0/hpoFreq[k] for k in hpoFreq}
    
    restartProb = 0.1
    
    BG_CALC = True
    GRAPH_BASED_WEIGHTS = True
    
    if BG_CALC:
        #background calculation
        bgNodes = mg.nodes['HPO']
        if GRAPH_BASED_WEIGHTS:
            #this equally distributes the weight across all terms
            bgProbs = {('HPO', h) : hpoWeights[h] for h in mg.nodes['HPO']}
            #bgProbs = {('HPO', h) : hpoWeights.get(h, 1.0) for h in mg.nodes['HPO']}
        else:
            #this distributes by how often we see weights
            bgProbs = {}
            '''
            for hpoTerms, geneSearch, codiRank, monarchRank, caseLabel in cases:
                for h in hpoTerms:
                    bgProbs[('HPO', h)] = bgProbs.get(('HPO', h), 0)+1#hpoWeights[h]
            '''
            '''
            fp = open('/Users/matt/data/UDN_HPO/MH_curated_HPO.csv', 'r')
            for l in fp:
                pieces = l.rstrip().split(',')
                if len(pieces) == 2:
                    ht = pieces[1].split(';')
                    for h in ht:
                        bgProbs[('HPO', h)] = bgProbs.get(('HPO', h), 0)+1.0/len(ht)#1.0#hpoWeights[h]
            fp.close()
            '''
            bgProbs = getHPODiseaseFreq()
            
            #for k in bgProbs.keys():
            #    bgProbs[k] *= hpoWeights[k[1]]
            
            #bgProbs = getPushHPODiseaseFreq('/Users/matt/data/HPO_dl/hp.obo')
        
        bg = mg.calculateBackground(bgProbs, restartProb)
    else:
        bg = None
    
    rankTypes = set(['gene'])
    
    for hpoTerms, geneSearch, codiRank, monarchRank, caseLabel in cases:
        #print caseLabel, hpoTerms
        if GRAPH_BASED_WEIGHTS:
            startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
            #startProbs = {('HPO', h) : hpoWeights.get(h, 1.0) for h in hpoTerms}
        else:
            startProbs = {('HPO', h) : 1.0 for h in hpoTerms}
            '''
            #this removes terms based on how much weight they distribute to their top 20 genes
            startProbs = {}
            for h in hpoTerms:
                tempRank = mg.RWR_rank({('HPO', h):1.0}, restartProb, rankTypes, bg)
                tw = sum([w for w, t, h in tempRank[0:20]])
                if tw > 0.1:
                    startProbs[('HPO', h)] = 1.0
            
            if len(startProbs) == 0:
                startProbs = {('HPO', h) : 1.0 for h in hpoTerms}
            '''
        #print startProbs
        rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)
        
        gs = set(codiRank)
        graphWeights = []
        for w, l, g in rankedGenes:
            if g in gs:
                graphWeights.append((w, g))
        
        for gene in sorted(geneSearch):
            if gene in codiRank:
                codiScore = codiRank.index(gene)
                monarchScore = (len(codiRank)-1+len(monarchRank)) // 2
                for x in range(0, len(monarchRank)):
                    if monarchRank[x][1] == gene:
                        w = monarchRank[x][0]
                        r = x
                        
                        #expand in the event of equals
                        l = r
                        h = r
                        while l > 0 and monarchRank[l-1][0] == w:
                            l -= 1
                        while h < len(monarchRank)-1 and monarchRank[h+1][0] == w:
                            h += 1
                        monarchScore = (l+h) // 2
                        break
                
                graphScore = (len(codiRank)-1+len(graphWeights)) // 2
                for x in range(0, len(graphWeights)):
                    if graphWeights[x][1] == gene:
                        w = graphWeights[x][0]
                        r = x
                        
                        #expand in the event of equals
                        l = r
                        h = r
                        while l > 0 and graphWeights[l-1][0] == w:
                            l -= 1
                        while h < len(graphWeights)-1 and graphWeights[h+1][0] == w:
                            h += 1
                        graphScore = (l+h) // 2
                        break
                
                print('\t'.join(str(v) for v in [caseLabel, ';'.join(hpoTerms), gene, codiScore, monarchScore, graphScore]))
            else:
                #for some reason, this gene isn't in the codi ranking
                pass
    
def rankCasesFromMana():
    retCases = []
    
    #load the genes for each case
    reportedGenes = '/Users/matt/data/codi_monarch_scores/UDN_reported_results.csv'
    caseToGenes = {}
    fp = open(reportedGenes, 'r')
    fp.readline()
    for l in fp:
        pieces = l.strip('\n').strip('\r').split(',')
        geneLabel = pieces[0]
        caseLabel = pieces[1]
        classification = pieces[11]
        
        if classification == 'primary':
            if caseLabel in caseToGenes:
                caseToGenes[caseLabel].add(geneLabel)
            else:
                caseToGenes[caseLabel] = set([geneLabel])
    fp.close()
    
    caseDirectories = sorted(glob.glob('/Users/matt/data/codi_monarch_scores/*/'))
    for cd in caseDirectories:
        caseLabel = cd.split('/')[-2]
        if len(caseToGenes.get(caseLabel, [])) == 0:
            print('no genes returned for '+caseLabel)
            continue
        
        #get all HPO terms for the case
        fn = cd+'HPO_input.txt'
        fp = open(fn, 'r')
        hpo = set([])
        for l in fp:
            if l[0:3] == 'HP:':
                hpo.add(l.rstrip())
        fp.close()
        
        #get the CODI rankings
        fn = cd+'codi_genes_parsed.csv'
        fp = open(fn, 'r')
        fp.readline()
        codiOrder = []
        for l in fp:
            pieces = l.rstrip().split(',')
            r, g, dummy = pieces
            codiOrder.append(g)
        
        #get the monarch rankings
        fn = cd+'out_monarch_scores.csv'
        fp = open(fn, 'r')
        fp.readline()
        monarchOrder = []
        for l in fp:
            pieces = l.rstrip().split(',')
            g, s, r = pieces
            monarchOrder.append((int(s), g))
        
        retCases.append((hpo, caseToGenes[caseLabel], codiOrder, monarchOrder, caseLabel))
    
    return retCases

def rankCase(mg):
    #training case 882
    #arthrogryposis, failure to thrive, pneumonia, retinopathy, ophthalmoplegia, dull facial expression, duane anomaly, abnormality of muscle size
    #hpoTerms = set(['HP:0002804', 'HP:0001508', 'HP:0002090', 'HP:0000580', 'HP:0000597', 'HP:0000338', 'HP:0009921', 'HP:0030236'])
    #minus things explained by piezo2
    #failure to thrive
    #hpoTerms = set(['HP:0001508'])
    
    #udn case 1
    #hpoTerms = set(['HP:0010286', 'HP:0002088', 'HP:0002643', 'HP:0000164', 'HP:0002110', 'HP:0000670', 'HP:0000217', 'HP:0000271', 'HP:0006536', 'HP:0006522', 'HP:0100755', 'HP:0000633', 'HP:0100759', 'HP:0002097'])
    #jsonDump = '/Users/matt/Downloads/SL126151_results.json'
    
    #udn case 2
    #hpoTerms = set(['HP:0000846', 'HP:0000252', 'HP:0000826', 'HP:0002444', 'HP:0000824', 'HP:0002902', 'HP:0001250', 'HP:0007165', 'HP:0000871'])
    #jsonDump = '/Users/matt/Downloads/SL126745_results.json'
    #jsonDump = '/Users/matt/Downloads/SL126745_results_f05.json'
    
    #udn case 3
    #hpoTerms = {'HP:0001510', 'HP:0002900', 'HP:0002164', 'HP:0002902', 'HP:0001802', 'HP:0001857', 'HP:0011968', 'HP:0005707'}
    #jsonDump = '/Users/matt/Downloads/SL126748_results.json'
    
    #udn case 15
    #hpoTerms = set(['HP:0002307', 'HP:0002141', 'HP:0002046', 'HP:0007286', 'HP:0001621', 'HP:0000722', 'HP:0000649', 'HP:0001268', 'HP:0012171', 'HP:0001152', 'HP:0002371', 'HP:0002376', 'HP:0002355', 'HP:0007089', 'HP:0002193', 'HP:0000020', 'HP:0002607', 'HP:0010845', 'HP:0006895', 'HP:0000651', 'HP:0200085', 'HP:0002015', 'HP:0000713', 'HP:0001290', 'HP:0011448', 'HP:0025121', 'HP:0001272', 'HP:0000565', 'HP:0001347'])
    #jsonDump = '/Users/matt/Downloads/SL139131_results.json'
    
    #udn case 24
    #hpoTerms = set(['HP:0000501', 'HP:0004322', 'HP:0000511', 'HP:0000648', 'HP:0200055', 'HP:0012385', 'HP:0001315', 'HP:0010498', 'HP:0001264', 'HP:0001263', 'HP:0012157', 'HP:0007340', 'HP:0003701', 'HP:0000252', 'HP:0000316', 'HP:0000358', 'HP:0007936', 'HP:0001773', 'HP:0003693', 'HP:0001290', 'HP:0000826', 'HP:0001007', 'HP:0005616', 'HP:0000486', 'HP:0001385', 'HP:0002650'])
    #jsonDump = '/Users/matt/Downloads/SL143666_results.json'
    
    #udn case 36
    #hpoTerms = set(['HP:0002716', 'HP:0030880', 'HP:0000520', 'HP:0001297', 'HP:0002729', 'HP:0003237', 'HP:0005978', 'HP:0100546', 'HP:0011109', 'HP:0000577', 'HP:0000821', 'HP:0002099', 'HP:0011834'])
    #jsonDump = '/Users/matt/Downloads/SL149780_results.json'
    
    #udn case 43
    #cleft palate, congenital hip dysplasia, muscle weakness, increase muscle fatigue, scoliosis, club feet, hypotonia, hypertonia, short stature, muscle hypoplasia, triangular face, brachycephaly, cupped ears, webbed neck, axillary pterygia, pectus excavatum, narrow chest, hypoplastic labia majora
    #hpoTerms = set(['HP:0000175', 'HP:0001385', 'HP:0001324', 'HP:0003750', 'HP:0002650', 'HP:0001762', 'HP:0001290', 'HP:0001276', 'HP:0004322', 'HP:0009004', 'HP:0000325', 'HP:0000248', 'HP:0000378', 'HP:0000465', 'HP:0001060', 'HP:0000767', 'HP:0000774', 'HP:0000059'])
    #jsonDump = '/Users/matt/Downloads/probando_uno_results.json'
    #jsonDump = '/Users/matt/Downloads/SL154670_results.json'
    
    #udn case 44
    #hpoTerms = set(['HP:0006610', 'HP:0000028', 'HP:0030084', 'HP:0010878', 'HP:0000954', 'HP:0000311', 'HP:0001537', 'HP:0001263', 'HP:0001265', 'HP:0002786', 'HP:0011451', 'HP:0000363', 'HP:0001276', 'HP:0000286', 'HP:0002421', 'HP:0002870', 'HP:0010720', 'HP:0100704', 'HP:0008872', 'HP:0000340', 'HP:0002360'])
    #jsonDump = '/Users/matt/Downloads/SL154671_results.json'
    
    #udn case 49
    #hpoTerms = set(['HP:0003691', 'HP:0000020', 'HP:0002607', 'HP:0011098', 'HP:0001250', 'HP:0000717', 'HP:0003323'])
    #jsonDump = '/Users/matt/Downloads/SL156321_results.json'
    
    #udn case 64
    #hpoTerms = set(['HP:0000548', 'HP:0003234', 'HP:0001508', 'HP:0001272', 'HP:0003124', 'HP:0002079', 'HP:0001320', 'HP:0001263', 'HP:0001251'])
    #jsonDump = '/Users/matt/Downloads/SL201104_results.json'
    
    #udn case 82
    #hpoTerms = set(['HP:0006530', 'HP:0005567', 'HP:0000121', 'HP:0000707', 'HP:0001644', 'HP:0006685'])
    #jsonDump = '/Users/matt/Downloads/SL209918_results.json'
    
    #udn case 84; second set is with the new info
    #hpoTerms = set(['HP:0000750', 'HP:0002451', 'HP:0001288', 'HP:0002080', 'HP:0006895', 'HP:0012531', 'HP:0003457', 'HP:0002141', 'HP:0007340', 'HP:0006315', 'HP:0040083', 'HP:0009763'])
    #hpoTerms = set(['HP:0000750', 'HP:0002451', 'HP:0001288', 'HP:0001009', 'HP:0001360', 'HP:0002080', 'HP:0006895', 'HP:0001258', 'HP:0003457', 'HP:0012531', 'HP:0002141', 'HP:0007340', 'HP:0006315', 'HP:0040083', 'HP:0009763'])
    #jsonDump = '/Users/matt/Downloads/SL210501_results.json'
    
    #udn case 107
    #hpoTerms = set(['HP:0000975', 'HP:0001290', 'HP:0000252', 'HP:0002079', 'HP:0001324', 'HP:0001263', 'HP:0004712', 'HP:0001321', 'HP:0002028', 'HP:0001508', 'HP:0001250', 'HP:0002650', 'HP:0002020', 'HP:0003128', 'HP:0003198'])
    #jsonDump = '/Users/matt/Downloads/SL219796_results.json'
    
    #udn case XXX
    hpoTerms = set(['HP:0003198', 'HP:0003202'])
    
    
    #primary data
    hpoWeights = pickle.load(open('/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle', 'rb'))
    startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
    restartProb = 0.1
    
    #background calculation
    bgNodes = mg.nodes['HPO']
    bgProbs = {('HPO', h) : hpoWeights[h] for h in mg.nodes['HPO']}
    bg = mg.calculateBackground(bgProbs, restartProb)
    #bg = None
    
    rankTypes = set(['gene'])
    
    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)
    
    print('Raw gene list:')
    for i, (w, l, g) in enumerate(rankedGenes[0:20]):
    #for i, (w, l, g) in enumerate(rankedGenes):
        print(i, w, l, g)
        #if g == 'PYROXD1':
        #    break
    print()
    
    codiGenes = getCodiGenes(jsonDump)
    print('Codi gene list:')
    count = 0
    for i, (w, l, g) in enumerate(rankedGenes):
        if g in codiGenes:
            print(count, w, l, g)
            count += 1
            if count >= 20:
                break

def getCodiGenes(jsonDump):
    dat = json.load(open(jsonDump, 'r'))
    ret = set([])
    for r in dat:
        geneList = r['genes']
        for v in geneList:
            for v2 in v['info']:
                if v2['label'] == 'gene symbol':
                    ret.add(v2['data'])
    return ret

def rankCasesFromUDN(mg):
    fp = open('/Users/matt/Downloads/MH_curated_HPO - MH_curated_HPO.tsv', 'r')
    dr = csv.DictReader(fp, delimiter='\t')
    
    siteResults = {}
    sites = set([])
    
    TERMINAL_CASE = 'case178'
    
    caseDict = {}
    caseOrder = []
    for d in dr:
        caseDict[d['case_label']] = d
        caseOrder.append(d['case_label'])
        
        s = d['site']
        sites.add(s)
        siteResults[(s, 'total')] = siteResults.get((s, 'total'), 0)+1
        for k in ['ret_5', 'ret_4', 'ret_3']:
            if d[k] != '':
                siteResults[(s, k)] = siteResults.get((s, k), 0)+1
                break
        
        for k in ['secondary', 'SV']:
            siteResults[(s, k)] = siteResults.get((s, k), 0)+(d[k] != '')
        
        if d['case_label'] == TERMINAL_CASE:
            break
            
    fp.close()
    
    print('SITE', 'PATH', 'LIKELY_PATH', 'VUS', 'SECONDARIES', 'SV', 'TOTAL', sep='\t')
    for s in sorted(sites):
        vals = [s,
                siteResults.get((s, 'ret_5'), 0),
                siteResults.get((s, 'ret_4'), 0),
                siteResults.get((s, 'ret_3'), 0),
                siteResults.get((s, 'secondary'), 0),
                siteResults.get((s, 'SV'), 0),
                siteResults.get((s, 'total'), 0)]
        
        for x in range(1, 6):
            vals[x] = vals[x] / vals[6]
        
        print(*vals, sep='\t')
    
    #primary data
    hpoWeights = pickle.load(open('/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle', 'rb'))
    restartProb = 0.1
    
    #background calculation
    bgNodes = mg.nodes['HPO']
    
    #'''
    #assume all are equal in background; NOTE: weights are scaled
    bgProbs = {('HPO', h) : hpoWeights[h] for h in mg.nodes['HPO']}
    #'''
    '''
    #calculate bg probs from case freq; NOTE: weights are 1.0
    bgProbs = {}
    for case in caseOrder:
        hpoTerms = caseDict[case]['HPO_terms'].split(';')
        for h in hpoTerms:
            bgProbs[('HPO', h)] = bgProbs.get(('HPO', h), 0)+1.0
    '''
    
    #calculate bg probs from disease frequencies
    #bgProbs = getHPODiseaseFreq()
    
    
    bg = mg.calculateBackground(bgProbs, restartProb)
    
    rankTypes = set(['gene'])
    
    results = {'ret_5':([], []),
               'ret_4':([], []),
               'ret_3':([], [])}
    
    allLines = []
    level = []
    rankings = []
    targets = []
    allTerms = []
    
    for case in caseOrder:
        hpoTerms = caseDict[case]['HPO_terms'].split(';')
        if len(hpoTerms) == 0:
            continue
        
        #sortedProbs = sorted([(hpoWeights[h], h) for h in hpoTerms], reverse=True)
        #startProbs = {('HPO', h) : hpoWeights[h] for w, h in sortedProbs[0:5]}
        '''
        startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
        #startProbs = {('HPO', h) : 1.0 for h in hpoTerms}
        rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)
        
        for k in ['ret_5', 'ret_4', 'ret_3']:
            if caseDict[case][k] != '':
                geneSet = set(caseDict[case][k].split(';'))
                print(case, geneSet)
                #print(startProbs)
                #print(sortedProbs)
                for i, (w, t, g) in enumerate(rankedGenes):
                    if (g in geneSet):
                        results[k][0].append(len(hpoTerms))
                        results[k][1].append(i)
                        print('', g, i, sep='\t')
        '''
        
        #startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
        for k in ['ret_5']:#, 'ret_4', 'ret_3']:
            if caseDict[case][k] != '':
                geneSet = set(caseDict[case][k].split(';'))
                print(case, geneSet)
                
                for gene in geneSet:
                    rankedTerms = []
                    for h in hpoTerms:
                        startProbs = {('HPO', h) : hpoWeights[h]}
                        rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)
                        tw = sum([w for w, t, g in rankedGenes[0:10]])
                        for i, (w, t, g) in enumerate(rankedGenes):
                            if gene == g:
                                #rankedTerms.append((i, h, tw))
                                rankedTerms.append((w, h, tw))
                                break
                        else:
                            #rankedTerms.append((i+1, h, tw))
                            rankedTerms.append((0.0, h, tw))
                    
                    #rankedTerms.sort()
                    rankedTerms.sort(reverse=True)
                    #print(rankedTerms[0:3])
                    
                    startProbs = {}
                    ordered = []
                    for w, h, tw in rankedTerms:
                        startProbs[('HPO', h)] = hpoWeights[h]
                        rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)
                        for i, (w, t, g) in enumerate(rankedGenes):
                            if gene == g:
                                v = i
                                break
                        else:
                            v = i+1
                        
                        ordered.append(v)
                    
                    allLines.append(ordered)
                    level.append(k)
                    rankings.append(rankedTerms)
                    targets.append(gene)
                    allTerms.append(v)
                    print('\t', gene, v)
    
    allTerms = np.array(allTerms)
    
    print()
    print('Mean:', np.mean(allTerms))
    print('Median:', np.median(allTerms))
    print('#>=20:', np.sum(allTerms >= 20))
    print()
    #exit()
    
    import matplotlib
    import matplotlib.pyplot as plt
        
    '''
    for k in results.keys():
        print(k, np.mean(results[k][1]), np.median(results[k][1]))
        
        plt.figure()
        plt.scatter(results[k][0], results[k][1])
        plt.show()
    '''
    
    SHOW_RANK_PLOT = True
    SHOW_SPEC_PLOT = False
    
    if SHOW_RANK_PLOT:
        cDict = {'ret_5':'red',
                 'ret_4':'blue',
                 'ret_3':'green'}
        plt.figure()
        for i, l in enumerate(allLines):
            am = np.argmin(l)
            
            plt.plot(range(1, len(l)+1), l, color=cDict[level[i]], alpha=0.1)
            plt.scatter([am+1], [l[am]], color=cDict[level[i]], alpha=0.2)
            
            print(targets[i], level[i], l[am])
            for x in range(0, am+1):
                print(rankings[i][x][2], rankings[i][x][1], l[x])
            print('final', 'final', l[-1])
            print()
            
        plt.yscale('symlog')
        plt.xlim(1, 10)
        plt.ylim(0)
        plt.grid()
        plt.xlabel('Number of Top HPO Terms')
        plt.ylabel('Ranking (points are absolute minima)')
        plt.show()
    
    if SHOW_SPEC_PLOT:
        pickleP2gFN = '/Users/matt/data/HPO_dl/multiP2G_biogrid_pushup.pickle'
        p2g = pickle.load(open(pickleP2gFN, 'rb'))
        
        cDict = {'ret_5':'red',
                 'ret_4':'blue',
                 'ret_3':'green'}
        plt.figure()
        bps = []
        for i, l in enumerate(allLines):
            am = np.argmin(l)
            #plt.plot(range(1, len(l)+1), l, color=cDict[level[i]], alpha=0.1)
            #plt.scatter([am+1], [l[am]], color=cDict[level[i]], alpha=0.2)
            
            if l[am] < 100:
                '''
                while len(bps) < am+1:
                    bps.append([])
                
                for j, (v, h) in enumerate(rankings[i][0:am+1]):
                    #bps[j].append(hpoWeights[h])
                    bps[j].append(len(p2g.get(h, {})))
                '''
                
                #vals = [len(p2g.get(h, {})) for j, (v, h, tw) in enumerate(rankings[i])]#[0:am+1])]
                #vals = [tw for j, (v, h, tw) in enumerate(rankings[i][0:am+1])]
                vals = [tw for j, (v, h, tw) in enumerate(rankings[i])]
                
                for x in range(0, am+1):
                    print(rankings[i][x][2], rankings[i][x][1])
                print()
                
                #plt.plot(range(1, am+2), vals, color=cDict[level[i]], alpha=0.2)
                #plt.plot(range(1, len(vals)+1), vals, color=cDict[level[i]], alpha=0.2)
                
                plt.plot(range(1, am+2), vals[0:am+1], color='blue', alpha=.2)
                plt.plot(range(am+1, len(vals)+1), vals[am:], color='red', alpha=.2)
                
        #plt.boxplot(bps)
        
        plt.xlim(1, 10)
        plt.ylim(0)
        #plt.yscale('symlog')
        plt.grid()
        plt.xlabel('HPO Term')
        plt.ylabel('Distribution of specificity')
        plt.show()
    
def getHPODiseaseFreq():
    bgFreq = {}
    fp = open('/Users/matt/data/HPO_dl/phenotype_annotation.tab', 'r')
    for l in fp:
        pieces = l.rstrip().split('\t')
        if pieces[0] == 'OMIM':
            hpoTerm = pieces[4]
            bgFreq[('HPO', hpoTerm)] = bgFreq.get(('HPO', hpoTerm), 0)+1
    fp.close()
    return bgFreq

def getPushHPODiseaseFreq(graphStructureFN):
    #second, build a graph
    nodes, edges, altIDMap = loadGraphStructure(graphStructureFN)
    parentDict = deriveParents(nodes, edges)
    
    omimToHPO = {}
    fp = open('/Users/matt/data/HPO_dl/phenotype_annotation.tab', 'r')
    for l in fp:
        pieces = l.rstrip().split('\t')
        if pieces[0] == 'OMIM':
            omimID = pieces[1]
            hpoTerm = pieces[4]
            
            if (omimID not in omimToHPO):
                omimToHPO[omimID] = set([])
            omimToHPO[omimID].add(hpoTerm)
            omimToHPO[omimID] |= parentDict.get(hpoTerm, set([]))
    fp.close()
    
    bgFreq = {}
    for omimID in omimToHPO:
        for h in omimToHPO[omimID]:
            bgFreq[('HPO', h)] = bgFreq.get(('HPO', h), 0)+1
    
    return bgFreq
    
if __name__ == '__main__':
    
    
    #p2gW = getP2GWeights()
    #exit()
    
    #static files we will be using
    hpoPhenoToGenoFN = '/Users/matt/data/HPO_dl/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'
    graphStructureFN = '/Users/matt/data/HPO_dl/hp.obo'
    
    #d = getPushHPODiseaseFreq(graphStructureFN)
    #for k in sorted(d.keys()):
    #    print(k, d[k])
    #exit()
    
    #output files
    REGENERATE_GRAPH = False
    pickleGraphFN = '/Users/matt/data/HPO_dl/multigraph.pickle'
    REGENERATE_WEIGHTS = False
    pickleHPOWeightFN = '/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle'
    pickleP2gFN = '/Users/matt/data/HPO_dl/multiP2G_biogrid_pushup.pickle'
    
    #load or generate the graph
    if (not REGENERATE_GRAPH) and os.path.exists(pickleGraphFN):
        print('Loading from "'+pickleGraphFN+'"')
        mg = pickle.load(open(pickleGraphFN, 'rb'))
    else:
        print('Generating new LayeredGraph')
        includeDisease, includeCPDB, includeOmnipath = False, False, False
        mg = createNewLayeredGraph(hpoPhenoToGenoFN, graphStructureFN, includeDisease, includeCPDB, includeOmnipath)
        
        print('Pickling...')
        fp = open(pickleGraphFN, 'wb+')
        pickle.dump(mg, fp)
        fp.close()
    
    #load or generate the HPO weights
    if (not REGENERATE_WEIGHTS) and os.path.exists(pickleHPOWeightFN):
        print('Loading from "'+pickleHPOWeightFN+'"')
        hpoWeights = pickle.load(open(pickleHPOWeightFN, 'rb'))
        p2g = pickle.load(open(pickleP2gFN, 'rb'))
    else:
        print('Generate new HPO weights')
        hpoWeights, p2g = createHPOWeights(hpoPhenoToGenoFN, graphStructureFN)
        
        print('Pickling...')
        fp = open(pickleHPOWeightFN, 'wb+')
        pickle.dump(hpoWeights, fp)
        fp.close()
        
        fp = open(pickleP2gFN, 'wb+')
        pickle.dump(p2g, fp)
        fp.close()
    
    #now run our tests
    print(mg)
    '''
    initNodes = {('HPO', 'HP:0003002'): 1.0}
    restartProb = .1
    retSet = set(['HPO', 'gene'])
    print(mg.RWR_rank(initNodes, restartProb, retSet)[0:10])
    #'''
    #'''
    st = time.time()
    #runTests(mg)
    testRankings(mg)
    print(time.time()-st)
    #'''
    
    #this is a singleton case we don't know the answer to
    #rankCase(mg)
    
    #run through all the UDN cases you curated and generate stats
    #rankCasesFromUDN(mg)
    