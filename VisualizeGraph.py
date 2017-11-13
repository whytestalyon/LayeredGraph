#!/usr/bin/env python
'''
VisualizeGraph.py

This file contains helpful sub-routines for generating images from a Random Walk run.
'''

import math
import os
import struct

from LayeredGraph import LayeredGraph

def saveGraphImage(mg, outFN, rankings=None, minWeight=0.0, drawEdgeWeights=False, nodeTypes=None):
    '''
    This function generates a dot file for graphviz to visualize the graph
    @param mg - the LayeredGraph to generate an image from
    @param outFN - the location to save the output (.dot is expected)
    @param rankings - the full rankings of all nodes in the graph (default: None, do not color the graph and visualize the whole graph)
    @param minWeight - the minimum weight from the ranking required to show up in the image (default: 0.0)
    @param drawEdgeWeights - if True, weight values will be included on the edges (default: False)
    '''
    #if we have ranks, create a dictionary of the weights for lookup later
    if rankings != None:
        rDict = {}
        for w, t, v in rankings:
            rDict[(t, v)] = w
    
    #open the file for writing
    fp = open(outFN, 'w+')
    fp.write('digraph food {\n')
    n = mg.nodes
    
    if nodeTypes == None:
        nodeTypes = sorted(n.keys())
    
    #iterate through all nodes in the graph
    #for k in sorted(n.keys()):
    for k in nodeTypes:
        for v in sorted(n[k]):
            if rankings == None:    
                #if there are no rankings, then always write the node
                fp.write(k+'_'+v+';\n')
            else:
                #we have rankings, so only write the node if it has sufficient weight
                r = rDict[(k, v)]
                if r < minWeight:
                    continue
                
                #all weights are in the range [0, 1], so scale that up to RGB 255 scale
                fc = int(math.floor(r*255))
                rgb = (255, 255-fc, 255-fc)
                fcHash = '#'+bytes.hex(struct.pack('BBB',*rgb))
                
                #write the node and include the weight
                fp.write('{}_{} [label="{}_{} ({:.4f})" style=filled fillcolor="{}"];\n'.format(k, v.replace(':', '_'), k, v.replace(':', '_'), r, fcHash));
            
            #now go through the nodes again looking for edges
            #for k2 in sorted(n.keys()):
            for k2 in nodeTypes:
                for v2 in sorted(n[k2]):
                    #make sure this node has enough weight to show up
                    if rankings != None:
                        r2 = rDict[(k2, v2)]
                        if r2 < minWeight:
                            continue
                    
                    #if an edges exists, it has a weight > 0
                    w = mg.getEdge(k, v, k2, v2)
                    if w > 0.0:
                        wn = mg.getEdge(k, v, k2, v2, True)
                        if drawEdgeWeights:
                            #include the raw weight and the normalized weight
                            #TODO: option for one or both?
                            fp.write('{}_{} -> {}_{} [label="{}({:.2f})"];\n'.format(k, v, k2, v2, w, wn))
                        else:
                            #only include the edge itself
                            fp.write('{}_{} -> {}_{};\n'.format(k, v.replace(':', '_'), k2, v2.replace(':', '_')))
                            
        
    fp.write('}\n')
    fp.close()

def visualize_RWR(dotPrefix, imagePrefix, mg, startProbs, restartProb, bg=None, cycleLimit=1000, minWeight=0.0):
    '''
    Run RWR and generate a dot file for each iteration.  Requires graphviz to be installed to run "dot".
    @param dotPrefix - dot files will be saved to <dotPrefix>.<iteration>.dot
    @param imagePrefix - image files will be saved to <imagePrefix>.<iteration>.png
    @param mg - an instance of LayeredGraph
    @param startProbs - same as LayeredGraph.RWR_rank(..)
    @param restartProb - same as LayeredGraph.RWR_rank(..)
    @param bg - same as LayeredGraph.RWR_rank(..)
    @param cycleLimit - same as LayeredGraph.RWR_rank(..)
    @param minWeight - the minimum weight on a node to visualize it (default: 0.0)
    '''
    #first, generate the iterator
    rankTypes = set(mg.nodes.keys())
    rwr_iter = mg.RWR_iter(startProbs, restartProb, rankTypes, bg, cycleLimit)
    for x, rankings in enumerate(rwr_iter):
        dotFN = '.'.join([dotPrefix, str(x), 'dot'])
        pngFN = '.'.join([imagePrefix, str(x), 'png'])
        
        #create the dot file, then run dot to generate the image file
        saveGraphImage(mg, dotFN, rankings, minWeight=minWeight, drawEdgeWeights=False)
        os.system('dot -Tpng -o '+pngFN+' '+dotFN)
    