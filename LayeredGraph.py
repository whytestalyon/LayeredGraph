#!/usr/bin/env python
'''
LayeredGraph.py

This file contains the implementation of a graph that allows for multiple "Layers" or "types" of nodes in the graph.  This also
contains the implementation of Random Walk with Restart (RWR) that can be used to generate relationships in the graph.
'''

import numpy as np
import random
from scipy.sparse import csr_matrix, lil_matrix
import time

class LayeredGraph:
    def __init__(self):
        #list of nodes initialized to nothing
        self.nodes = {}
        
        #if noMultiWeights is True, then typeTransProb will do nothing
        self.uniformGraph = True
        self.typeTransProbs = {}
        
    def addNode(self, nodeType, nodeName):
        '''
        Adds a node to the graph
        @param nodeType - the node type (i.e. sub-graph label)
        @param nodeName - the label of the node in its subgraph
        '''
        if nodeType not in self.nodes:
            self.nodes[nodeType] = set([])
        self.nodes[nodeType].add(nodeName)
    
    def finalizeNodeList(self):
        '''
        This function finalizes the node list and constructs the edge matrix to be filled in by addEdge(...).  NOTE: Call this function only when
        ALL nodes that are required have been added using addNode(...).
        '''
        #all nodes have been added so initialize these values 
        self.nodeIndex = {}
        self.nodeTypeStart = {}
        self.nodeTypeLen = {}
        self.nodeOrdered = []
        
        #now calculate the index of each node in our matrix
        x = 0
        for nt in sorted(self.nodes.keys()):
            self.nodeTypeStart[nt] = x
            self.nodeTypeLen[nt] = len(self.nodes[nt])
            
            for nl in sorted(self.nodes[nt]):
                self.nodeIndex[(nt, nl)] = x
                self.nodeOrdered.append((nt, nl))
                x += 1
        
        #this will be indexed by [source, dest]
        self.edgeWeights = lil_matrix((x, x), dtype='float')
        
    def addEdge(self, st, sl, dt, dl, weight, isUndirected):
        '''
        This function adds an edge to the graph.
        @param st - source node type
        @param sl - source node label
        @param dt - destination node type
        @param dl - destination node label
        @param weight - the edge weight, use 1.0 if all edges are identical (i.e. unweighted)
        @param isUndirected - if True, both the described edge and its reverse will be added (source and dest are interchangeable)
        '''
        sourceInd = self.nodeIndex[(st, sl)]
        destInd = self.nodeIndex[(dt, dl)]
        
        #always add the directed edge, add the undirected if told to do so
        self.edgeWeights[sourceInd, destInd] = weight
        if isUndirected:
            self.edgeWeights[destInd, sourceInd] = weight
        
    def getEdge(self, st, sl, dt, dl, normalized=False):
        '''
        This function returns the weight on an edge in the graph
        @param st - source node type
        @param sl - source node label
        @param dt - destination node type
        @param dl - destination node label
        @param normalized - if True, return the normalized weight value [0.0, 1.0] from the computation; if False, return the raw weight from the addEdge(...) function
        @return - if the edge exists, the numerical weight value is returned; otherwise, this returns 0.0 to indicate an absent edge
        '''
        sourceInd = self.nodeIndex[(st, sl)]
        destInd = self.nodeIndex[(dt, dl)]
        if normalized:
            #this matrix was transposed, so flip the index order
            return self.transitionProbs[destInd, sourceInd]
        else:
            #raw weight values
            return self.edgeWeights[sourceInd, destInd]
    
    def setGraphJumpEqual(self):
        '''
        This function calculates an equal weight between graph jumping such that if there are N nodetypes, there is a 1/N probability of jumping
        to any other graph (including the current) as long as at least one edges exists that will enable that jump.  If an edge doesn't exist, it
        is handled in the ranking algorithm.
        '''
        self.uniformGraph = False
        for nt in list(self.nodes.keys()):
            for nt2 in list(self.nodes.keys()):
                self.typeTransProbs[(nt, nt2)] = 1.0/len(list(self.nodes.keys()))
    
    def calculateTransitionMatrix(self):
        '''
        This function will calculate the transition matrix based on all previously input data.  NOTE: Only call this after all edges have been added to
        the graph via addEdge(...).  
        '''
        #transition probs will be indexed by [dest, source] due to math reasons
        self.transitionProbs = self.edgeWeights.transpose()
        
        #1 - make sure every row has a non-zero sum; any that do automatically get a self-looping edge
        cs = np.sum(self.transitionProbs, axis=0)
        csz = np.where(cs == 0)[1]
        self.transitionProbs[csz, csz] = 1.0
        
        #this generally speeds everything up
        self.transitionProbs = csr_matrix(self.transitionProbs)
        
        #TODO: normalization takes the longest time of any steps by far.  Is there a smarter way to do it?
        #print 'normalize'
        #st = time.time()
        
        #2 - perform w/e normalization was specified
        if self.uniformGraph:
            self.normalizeUniformally()
        else:
            self.normalizeTransProbs()
        
        #print time.time()-st
        
    def normalizeUniformally(self):
        '''
        This function normalizes the transition probabilities to sum to 1.  This version treat all nodes and edges uniformally.
        There is no distinguishing between node types (although weight will impact the probabilities). 
        '''
        cs = np.sum(self.transitionProbs, axis=0)
        self.transitionProbs /= cs
    
    def normalizeTransProbs(self):
        '''
        This function normalizes the transition probabilities to sum to 1.  This version factors in graph jump probabilities stored in self.typeTransProbs.
        This approach allows us to control the rate at which probabilities move from one sub-graph type to another.  It requires self.typeTransProbs to have
        a value for EVERY (node type, node type) pair.
        '''
        #we have to handle each source (column) separately because we cannot guarantee edges always exist
        for ci in range(0, self.transitionProbs.shape[0]):
            st = self.nodeOrdered[ci][0]
            
            #go through and pull out transition probabilities if it's not empty
            tps = []
            tots = []
            for dt in sorted(self.nodes.keys()):
                dtStart = self.nodeTypeStart[dt]
                dtEnd = dtStart+self.nodeTypeLen[dt]
                s = np.sum(self.transitionProbs[dtStart:dtEnd, ci])
                tots.append(s)
                if s == 0:
                    tps.append(0.0)
                else:
                    tps.append(self.typeTransProbs[(st, dt)])
            
            #normalize based on what data is available
            tps = np.array(tps)
            tps /= np.sum(tps)
            
            for i, dt in enumerate(sorted(self.nodes.keys())):
                if tps[i] != 0:
                    dtStart = self.nodeTypeStart[dt]
                    dtEnd = dtStart+self.nodeTypeLen[dt]
                    self.transitionProbs[dtStart:dtEnd, ci] *= tps[i] / tots[i]
        
    def RWR_rank(self, startProbs, restartProb, rankTypes, bg=None, cycleLimit=1000):
        '''
        Random Walk with Restart method
        @param startProbs - the non-zero probabilities of restarting at a node
        @param restartProb - the probability of having a restart
        @param rankTypes - set of nodeTypes to be included in the return
        @param bg - background levels; if None, no effect; otherwise, must be a numpy array of the same size of the probabilities that is used to normalize before ranking
        @param cycleLimit - the maximum number of matrix multiplications to perform (default: 1000)
        @return - a list of tuples of form (weight, nodeType, nodeLabel) sorted with largest weights first
            weight - the probability of being at the node
            nodeType - type of the node, must be a member of "rankTypes" set
            nodeLabel - label of the node
        '''
        assert(restartProb >= 0.0 and restartProb <= 1.0)
        
        #create the initial probabilities (based on a restart)
        currentProbs = np.zeros(shape=(self.transitionProbs.shape[0], ), dtype='float')
        for k in list(startProbs.keys()):
            currentProbs[self.nodeIndex[k]] = startProbs[k]
        currentProbs /= np.sum(currentProbs)
        initProbs = np.copy(currentProbs)
        
        tm = self.transitionProbs
        
        assert(np.all(np.isclose(np.sum(tm, axis=0), 1)))
        assert(np.isclose(np.sum(currentProbs), 1))
        
        prev = np.zeros(shape=currentProbs.shape)
        x = 0
        while x < cycleLimit and not np.all(np.isclose(currentProbs, prev)):
            prev = currentProbs
            
            #classic formula: prob_(t+1) = (1-r)*TM.dot(prob_t) + r*prob_0
            currentProbs = (1-restartProb)*tm.dot(prev)+restartProb*initProbs
            x += 1
        
        #if we include background values, subtract them out to normalize the data
        if not (bg is None):
            currentProbs -= bg
        
        #now rank everything and only return values that are in the return set
        cpSort = np.argsort(currentProbs)[::-1]
        ret = []
        for v in cpSort:
            nt, nl = self.nodeOrdered[v]
            if nt in rankTypes:
                ret.append((currentProbs[v], nt, nl))
        
        return ret
    
    def calculateBackground(self, startProbs, restartProb, cycleLimit=1000):
        '''
        Random Walk with Restart background calculation.  This function performs a RWR using the given start and restart probabilities.  It then
        returns a blanket array of final probabilities for each node in the entire graph.  The intent is to provide a mechanism for normalizing
        the data by passing in start probabilities that have an equal likelihood of starting anywhere in the graph (or maybe a subgraph of known
        possible inputs).  This background can then be passed into the main RWR_rank method to calculate the relative increase/decrease in likelihood
        of arriving at a particular node when compared to the background.
        @param startProbs - the non-zero probabilities of restarting at a node, must be an entry for EVERY node that could be a source node
        @param restartProb - the probability of having a restart
        @param cycleLimit - the maximum number of matrix multiplications to perform (default: 1000)
        @return - an array containing the final probabilities for these conditions
        '''
        assert(restartProb >= 0.0 and restartProb <= 1.0)
        
        #create the initial probabilities (based on a restart)
        currentProbs = np.zeros(shape=(self.transitionProbs.shape[0], ), dtype='float')
        for k in list(startProbs.keys()):
            currentProbs[self.nodeIndex[k]] = startProbs[k]
        currentProbs /= np.sum(currentProbs)
        initProbs = np.copy(currentProbs)
        
        tm = self.transitionProbs
        
        assert(np.all(np.isclose(np.sum(tm, axis=0), 1)))
        assert(np.isclose(np.sum(currentProbs), 1))
        
        prev = np.zeros(shape=currentProbs.shape)
        x = 0
        while x < cycleLimit and not np.all(np.isclose(currentProbs, prev)):
            prev = currentProbs
            
            #classic formula: prob_(t+1) = (1-r)*TM.dot(prob_t) + r*prob_0
            currentProbs = (1-restartProb)*tm.dot(prev)+restartProb*initProbs
            x += 1
        
        return currentProbs
    
    def RWR_iter(self, startProbs, restartProb, rankTypes, bg=None, cycleLimit=1000):
        '''
        Random Walk with Restart method
        @param startProbs - the non-zero probabilities of restarting at a node
        @param restartProb - the probability of having a restart
        @param rankTypes - set of nodeTypes to be included in the return
        @param bg - background levels; if None, no effect; otherwise, must be a numpy array of the same size of the probabilities that is used to normalize before ranking
        @param cycleLimit - the maximum number of matrix multiplications to perform (default: 1000)
        @return - a list of tuples of form (weight, nodeType, nodeLabel) sorted with largest weights first
            weight - the probability of being at the node
            nodeType - type of the node, must be a member of "rankTypes" set
            nodeLabel - label of the node
        '''
        assert(restartProb >= 0.0 and restartProb <= 1.0)
        
        #create the initial probabilities (based on a restart)
        currentProbs = np.zeros(shape=(self.transitionProbs.shape[0], ), dtype='float')
        for k in list(startProbs.keys()):
            currentProbs[self.nodeIndex[k]] = startProbs[k]
        currentProbs /= np.sum(currentProbs)
        initProbs = np.copy(currentProbs)
        
        tm = self.transitionProbs
        
        assert(np.all(np.isclose(np.sum(tm, axis=0), 1)))
        assert(np.isclose(np.sum(currentProbs), 1))
        
        prev = np.zeros(shape=currentProbs.shape)
        x = 0
        while x < cycleLimit and not np.all(np.isclose(currentProbs, prev)):
            #now set the probs
            prev = np.copy(currentProbs)
            
            #if we include background values, subtract them out to normalize the data
            if not (bg is None):
                currentProbs -= bg
            
            #now rank everything and only return values that are in the return set
            cpSort = np.argsort(currentProbs)[::-1]
            ret = []
            for v in cpSort:
                nt, nl = self.nodeOrdered[v]
                if nt in rankTypes:
                    ret.append((currentProbs[v], nt, nl))
            yield ret
            
            #classic formula: prob_(t+1) = (1-r)*TM.dot(prob_t) + r*prob_0
            currentProbs = (1-restartProb)*tm.dot(prev)+restartProb*initProbs
            x += 1
        
        #if we include background values, subtract them out to normalize the data
        if not (bg is None):
            currentProbs -= bg
        
        #now rank everything and only return values that are in the return set
        cpSort = np.argsort(currentProbs)[::-1]
        ret = []
        for v in cpSort:
            nt, nl = self.nodeOrdered[v]
            if nt in rankTypes:
                ret.append((currentProbs[v], nt, nl))
        
        yield ret
    
if __name__ == '__main__':
    #example graph constructed and run here
    mg = LayeredGraph()
    
    mg.addNode('t1', 'n1')
    mg.addNode('t1', 'n2')
    mg.addNode('t2', 'n1')
    mg.addNode('t2', 'n2')
    mg.addNode('t3', 'ABA')
    mg.addNode('t3', 'C1orf72')
    
    mg.finalizeNodeList()
    
    mg.addEdge('t1', 'n1', 't2', 'n1', 1, True)
    mg.addEdge('t1', 'n1', 't2', 'n2', 1, True)
    mg.addEdge('t1', 'n2', 't2', 'n1', 1, True)
    #mg.addEdge('t1', 'n2', 't2', 'n2', 1, True)
    mg.addEdge('t2', 'n1', 't2', 'n2', 1, True)
    mg.addEdge('t2', 'n1', 't3', 'ABA', 1, False)
    mg.addEdge('t2', 'n2', 't3', 'C1orf72', 1, False)
    
    mg.setGraphJumpEqual()
    mg.calculateTransitionMatrix()
    
    ranks = mg.RWR_rank({('t1', 'n1'): 1.0, ('t1', 'n2'): 1.0}, 0.1, set(['t2', 't3']))
    print(ranks)
    
    print(mg.getEdge('t1', 'n1', 't2', 'n1'), mg.getEdge('t1', 'n1', 't2', 'n1', True))
    print(mg.getEdge('t1', 'n1', 't1', 'n1'), mg.getEdge('t1', 'n1', 't1', 'n1', True))
    
    