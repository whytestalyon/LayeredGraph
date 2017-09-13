'''
This will server as a basic wrapper for the LayeredGraph.  It will primarily be an API with minimal actual HTML utility.
'''

from flask import Flask
from flask import g
from flask import jsonify
from flask import render_template
from flask import Response
from flask import request
import pickle

from HPOParser import getTermsAndSyns

app = Flask(__name__)

@app.route('/')
def hello_world():
    ret = 'Hello, World!'
    ret += '\n'+str(mg)
    return ret

@app.route('/background')
def background():
    '''
    TODO: Function to return the background array, can be removed later OR folded into the initialization so background isn't constantly recalculated
    '''
    #primary data
    mg = getMultigraphVars()
    hpoWeights = pickle.load(open('/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle', 'r'))
    restartProb = 0.1
    
    #background calculation
    bgNodes = mg.nodes['HPO']
    bgProbs = {('HPO', h) : hpoWeights[h] for h in mg.nodes['HPO']}
    bg = mg.calculateBackground(bgProbs, restartProb)
    
    return jsonify(bg)

@app.route('/search')
def search():
    terms = getTermsAndSyns('/Users/matt/data/HPO_dl/hp.obo')
    return render_template('search.html', terms=terms)

@app.route('/rank', methods=['POST'])
def rank():
    '''
    This is intended to be the workhorse function.  POST needs to include a "term" list that should be HPO terms matching values from the layered graph.
    '''
    #"constant" global data
    mg, restartProb, hpoWeights, bg = getMultigraphVars()
    
    #TODO: make this actually work
    hpoTerms = set(request.form.getlist('term[]'))
    #startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
    startProbs = {}
    usedTerms = set([])
    missingTerms = set([])
    for h in hpoTerms:
        if hpoWeights.has_key(h):
            startProbs[('HPO', h)] = hpoWeights[h]
            usedTerms.add(h)
        else:
            missingTerms.add(h)
    
    rankTypes = set(['gene'])
    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)
    
    rankings = []
    for w, t, l in rankedGenes:
        rankings.append({'weight':w, 'nodeType':t, 'label':l})
    
    ret = {'rankings': rankings,
           'usedTerms': list(usedTerms),
           'missingTerms': list(missingTerms)}
    
    return jsonify(ret)
    
def getMultigraphVars():
    '''
    Get any variables that are constant between requests from flask.g.
    Note: I don't think this is doing quite what I want it to because the files are being loaded with each request, but it is relatively fast so moving on for now.
    @return 
    '''
    if not hasattr(g, 'mg'):
        #load or generate the graph
        pickleGraphFN = '/Users/matt/data/HPO_dl/multigraph.pickle'
        print 'Loading from "'+pickleGraphFN+'"'
        
        #these are variable that need to be set for multiple uses
        g.mg = pickle.load(open(pickleGraphFN, 'r'))
    
    if not hasattr(g, 'restartProb'):
        g.restartProb = 0.1
    
    if not hasattr(g, 'hpoWeights'):
        g.hpoWeights = pickle.load(open('/Users/matt/data/HPO_dl/multiHpoWeight_biogrid_pushup.pickle', 'r'))
    
    if not hasattr(g, 'bg'):
        bgProbs = {('HPO', h) : g.hpoWeights[h] for h in g.mg.nodes['HPO']}
        g.bg = g.mg.calculateBackground(bgProbs, g.restartProb)
    
    return (g.mg, g.restartProb, g.hpoWeights, g.bg)
    
# run the app.
if __name__ == "__main__":
    # Setting debug to True enables debug output. This line should be
    # removed before deploying a production app.
    
    #initializeLayeredGraph()
    
    app.debug = True
    app.run()