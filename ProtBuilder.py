#!/usr/bin/env python3
'''
Usage: python3 ProtBuilder.py

This file contains the methods to construct the LayeredGraph used by PyxisMap.  Pathing is currently local to Matt Holt's desktop.
'''

import json
import numpy as np
import pickle

import LayeredGraph

def createProtGraph():
    '''
    This function will create an the graph we use for protein protein interactions
    @return - an instance of LayeredGraph containing the constructed graph
    '''
    #first, load the dictionary that pairs proteins with genes
    prot2gene = loadProteinDict()
    print('Loaded '+str(len(prot2gene))+' protein-gene pairs.')
    
    #second, load the protein-protein interactions
    ppi = loadPPI()
    print('Loaded '+str(len(ppi))+' protein-protein interaction pairs.')
    
    #create the graph
    mg = LayeredGraph.LayeredGraph()
    
    #add all genes that have an entry in the prot-prot interactions
    print('Adding gene nodes...')
    for p1, p2, conf, participants in ppi:
        for prot in p1.split('.'):
            mg.addNode('gene', prot2gene.get(prot, prot))
        for prot in p2.split('.'):
            mg.addNode('gene', prot2gene.get(prot, prot))
        
    #finalize graph nodes
    mg.finalizeNodeList()
    
    #add all the prot-prot interactions as edges
    print('Calculating edge weights...')
    edgeWeights = {}
    for p1, p2, conf, participants in ppi:
        #TODO: change this to be a function of confidence & # of participants
        if conf != 'NA':
            #interactionWeight = 1.0
            interactionWeight = conf
            
            for prot1 in p1.split('.'):
                g1 = prot2gene.get(prot1, prot1)
                for prot2 in p2.split('.'):
                    g2 = prot2gene.get(prot2, prot2)
                    
                    #add in the sorted order
                    if g1 < g2:
                        edgeWeights[(g1, g2)] = max(interactionWeight, edgeWeights.get((g1, g2), 0.0))
                    else:
                        edgeWeights[(g2, g1)] = max(interactionWeight, edgeWeights.get((g2, g1), 0.0))
    
    #add all the weights to the graph now
    print('Adding gene-gene edges based on protein-protein interactions...')
    for g1, g2 in edgeWeights.keys():
        mg.addEdge('gene', g1, 'gene', g2, edgeWeights[(g1, g2)], True)
    
    #these are the final steps
    print('Calculating final transition matrix...')
    mg.setGraphJumpEqual()
    mg.calculateTransitionMatrix()
    
    return mg

def loadProteinDict():
    '''
    This function will read in the file that pairs protein labels (from Uniprot?) with gene names
    @return - a dictionary where key is the protein label and value is the gene name
    '''
    fn = '/Users/matt/data/HPO_dl/uniprot-all.tab'
    fp = open(fn, 'r')
    fp.readline()
    
    ret = {}
    for l in fp:
        pieces = l.rstrip().split('\t')
        
        uniprotLabel = pieces[1]
        status = pieces[2]
        geneNames = pieces[4]
        organism = pieces[5]
        
        assert(status == 'reviewed')
        assert(organism == 'Homo sapiens (Human)')
        assert(uniprotLabel not in ret)
        
        ret[uniprotLabel] = geneNames.split(' ')[0]
    
    fp.close()
    return ret

def loadPPI():
    '''
    This function will read in the file that contains protein-protein interactions and return the information we care about
    @return - a list of tuples where each tuple contains (protein1, protein2, confidence, number of participants) such that protein1 lexicographically precedes protein2
    '''
    fn = '/Users/matt/data/HPO_dl/ConsensusPathDB_human_PPI'
    fp = open(fn, 'r')
    fp.readline()
    fp.readline()
    
    ret = []
    for l in fp:
        pieces = l.rstrip().split('\t')
        participants = sorted(pieces[2].split(','))
        confidence = pieces[3]
        
        if len(participants) == 1:
            #currently, we do not care about these
            pass
        else:
            confidence = (float(confidence) if confidence != 'NA' else 'NA')
            for x in range(0, len(participants)):
                for y in range(x+1, len(participants)):
                    ret.append((participants[x], participants[y], confidence, len(participants)))
    
    fp.close()
    return ret

if __name__ == '__main__':
    #TODO: maybe write code/scripts to fetch these? I know the local HPO datasets are fairly old now
    
    #static files we will be using
    
    #output files
    pickleGraphFN = '/Users/matt/data/HPO_dl/protgraph.pickle'
    
    #create the graph
    mg = createProtGraph()
    
    #load or generate the graph
    print('Pickling ProtGraph...')
    fp = open(pickleGraphFN, 'wb+')
    pickle.dump(mg, fp)
    fp.close()
    
    