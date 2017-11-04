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
import requests
from urllib.parse import urlencode, quote_plus
from pprint import pprint
from HPOParser import getTermsAndSyns

app = Flask(__name__)

@app.route('/background')
def background():
    '''
    TODO: Function to return the background array, can be removed later OR folded into the initialization so background isn't constantly recalculated
    '''
    #primary data
    mg = getMultigraphVars()
    hpoWeights = pickle.load(open('./HPO_graph_data/multiHpoWeight_biogrid_pushup.pickle', 'rb'))
    restartProb = 0.1
    
    #background calculation
    bgNodes = mg.nodes['HPO']
    bgProbs = {('HPO', h) : hpoWeights[h] for h in mg.nodes['HPO']}
    bg = mg.calculateBackground(bgProbs, restartProb)
    
    return jsonify(bg)

@app.route('/')
@app.route('/search')
def search():
    terms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    return render_template('search.html', terms=terms)

@app.route('/text')
def text():
    return render_template('textToRank.html')

@app.route('/rank', methods=['POST'])
def rank():
    '''
    This is intended to be the workhorse function.  POST needs to include a "term" list that should be HPO terms matching values from the layered graph.
    '''
    # "constant" global data
    mg, restartProb, hpoWeights, bg = getMultigraphVars()

    # TODO: make this actually work
    hpoTerms = set([str(x) for x in request.form.getlist('term[]')])
    pprint(hpoTerms)
    # startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
    startProbs = {}
    usedTerms = set([])
    missingTerms = set([])
    for h in hpoTerms:
        if h in hpoWeights:
            startProbs[('HPO', h)] = hpoWeights[h]
            usedTerms.add(h)
        else:
            missingTerms.add(h)

    rankTypes = set(['gene'])
    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)

    if request.form['action'] == 'JSON':
        rankings = []
        for w, t, l in rankedGenes:
            rankings.append({'weight': w, 'nodeType': t, 'label': l})

        ret = {'rankings': rankings,
               'usedTerms': list(usedTerms),
               'missingTerms': list(missingTerms)}

        return jsonify(ret)
    else:
        ret = 'Used: ' + str(usedTerms) + '<br>'
        ret += 'CSV-used: ' + str(';'.join(sorted(usedTerms))) + '<br>'
        ret += 'Missing: ' + str(missingTerms) + '<br><br>'
        for i, (w, t, l) in enumerate(rankedGenes[0:20]):
            ret += ' '.join([str(x) for x in (i, w, t, l)]) + '<br>'
        return ret@app.route('/rank', methods=['POST'])

@app.route('/text/rank', methods=['POST'])
def textrank():
    '''
    Workhorse function for converting text to hpo terms to genes.  POST needs to include a "indications" payload that
    should be text to be annotated by the NCBO annotator.
    '''

    indicationsText = request.form['indications']
    payload = {'apikey': '8b5b7825-538d-40e0-9e9e-5ab9274a9aeb',
               'text': indicationsText,
               'ontologies': 'HP',
               'whole_word_only': 'false'}
    url = 'http://data.bioontology.org/annotator?' + urlencode(payload, quote_via=quote_plus)
    resp = requests.get(url)
    pprint(resp.json())

    # "constant" global data
    mg, restartProb, hpoWeights, bg = getMultigraphVars()

    # TODO: make this actually work
    hpoTerms = set([])
    for annot in resp.json():
        hpoTerms.add(annot['annotatedClass']['@id'].split('/')[-1].replace('_', ':'))
        # hpoTerms.add(annot['annotatedClass']['@id'].replace('_', ':'))

    pprint(hpoTerms)
    # startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
    startProbs = {}
    usedTerms = set([])
    missingTerms = set([])
    for h in hpoTerms:
        if h in hpoWeights:
            startProbs[('HPO', h)] = hpoWeights[h]
            usedTerms.add(h)
        else:
            missingTerms.add(h)

    rankTypes = set(['gene'])
    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)

    if request.form['action'] == 'JSON':
        rankings = []
        for w, t, l in rankedGenes:
            rankings.append({'weight': w, 'nodeType': t, 'label': l})

        ret = {'rankings': rankings,
               'usedTerms': list(usedTerms),
               'missingTerms': list(missingTerms)}

        return jsonify(ret)
    elif request.form['action'] == 'GENES':
        ret = ''
        for w, t, l in rankedGenes:
            ret += l + ';'
        return ret
    else:
        ret = 'Used: ' + str(usedTerms) + '<br>'
        ret += 'CSV-used: ' + str(';'.join(sorted(usedTerms))) + '<br>'
        ret += 'Missing: ' + str(missingTerms) + '<br><br>'
        for i, (w, t, l) in enumerate(rankedGenes[0:20]):
            ret += ' '.join([str(x) for x in (i, w, t, l)]) + '<br>'
        return ret

def getMultigraphVars():
    '''
    Get any variables that are constant between requests from flask.g.
    Note: I don't think this is doing quite what I want it to because the files are being loaded with each request, but it is relatively fast so moving on for now.
    @return 
    '''
    if not hasattr(g, 'mg'):
        #load or generate the graph
        pickleGraphFN = './HPO_graph_data/multigraph.pickle'
        print('Loading from "'+pickleGraphFN+'"')
        
        #these are variable that need to be set for multiple uses
        g.mg = pickle.load(open(pickleGraphFN, 'rb'))
    
    if not hasattr(g, 'restartProb'):
        g.restartProb = 0.1
    
    if not hasattr(g, 'hpoWeights'):
        g.hpoWeights = pickle.load(open('./HPO_graph_data/multiHpoWeight_biogrid_pushup.pickle', 'rb'))
    
    if not hasattr(g, 'bg'):
        bgProbs = {('HPO', h) : g.hpoWeights[h] for h in g.mg.nodes['HPO']}
        g.bg = g.mg.calculateBackground(bgProbs, g.restartProb)
    
    return (g.mg, g.restartProb, g.hpoWeights, g.bg)
    
# run the app.
if __name__ == "__main__":
    # Setting debug to True enables debug output. This line should be
    # removed before deploying a production app.
    
    #initializeLayeredGraph()

    app.run(debug=True, host='0.0.0.0')