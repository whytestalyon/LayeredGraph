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
import xmltodict
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

@app.route('/search')
def search():
    '''
    DEPRECATED
    This is the original search page MH used
    '''
    terms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    return render_template('search.html', terms=terms)

@app.route('/')
@app.route('/text', methods=['GET'])
def text():
    '''
    This is the main page that a user will hit that only performs a search over all terms
    '''
    return render_template('textToRank.html')

@app.route('/table', methods=['GET'])
def table():
    '''
    This is the page for the table view, it performs a more complex query so it is given its own page
    '''
    return render_template('textToTable.html')

@app.route('/terms', methods=['GET'])
def terms():
    '''
    This function returns a JSON presentation of the terms that match user typed text
    '''
    pprint(request.args)
    searchTerm = str(request.args.get('term'))
    hpoTerms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    options = list([])

    for hpo, defText in hpoTerms:
        if searchTerm.lower() in defText.lower():
            options.append({'id': hpo, 'text': (hpo + ' ' + defText)})
        if searchTerm.lower() in hpo.lower():
            options.append({'id': hpo, 'text': (hpo + ' ' + defText)})

    return jsonify({'results': options})

@app.route('/rank', methods=['POST'])
def rank():
    '''
    This is intended to be the workhorse function.  POST needs to include a "term" list that should be HPO terms matching values from the layered graph.
    '''
    mydata = request.get_json()

    # "constant" global data
    mg, restartProb, hpoWeights, bg = getMultigraphVars()
    hpoTerms = set([str(x) for x in mydata])

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

    rankings = []
    for w, t, l in rankedGenes:
        rankings.append({'weight': w, 'nodeType': t, 'label': l})

    ret = {'rankings': rankings,
           'usedTerms': list(usedTerms),
           'missingTerms': list(missingTerms)}
    return jsonify(ret)

@app.route('/deeprank', methods=['POST'])
def deeprank():
    '''
    This is intended to be the workhorse function.  POST needs to include a "term" list that should be HPO terms matching values from the layered graph.
    '''
    mydata = request.get_json()

    # "constant" global data
    mg, restartProb, hpoWeights, bg = getMultigraphVars()
    hpoTerms = set([str(x) for x in mydata])
    rankTypes = set(['gene'])
    
    # startProbs = {('HPO', h) : hpoWeights[h] for h in hpoTerms}
    startProbs = {}
    usedTerms = set([])
    missingTerms = set([])
    
    indivRanks = {}
    
    for h in hpoTerms:
        if h in hpoWeights:
            startProbs[('HPO', h)] = hpoWeights[h]
            usedTerms.add(h)
            
            #do an individual weighting for each term as we go
            termWeights = mg.RWR_rank({('HPO', h) : 1.0}, restartProb, rankTypes, bg)
            for i, (w, t, l) in enumerate(termWeights):
                if (l not in indivRanks):
                    indivRanks[l] = {}
                indivRanks[l][h] = (i+1, w)
                #residual method for when these were separate values
                #indivRanks[l][h] = w
                #indivRanks[l][h+'-RANK'] = (i+1)
        else:
            missingTerms.add(h)

    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)

    rankings = []
    for i, (w, t, l) in enumerate(rankedGenes):
        temp = {'weight': w, 'nodeType': t, 'label': l, 'rank': i+1}
        temp.update(indivRanks[l])
        rankings.append(temp)

    ret = {'rankings': rankings,
           'usedTerms': list(usedTerms),
           'missingTerms': list(missingTerms)}
    return jsonify(ret)

@app.route('/text/annotate', methods=['GET'])
def textannotate():
    '''
    This function hits the NIH annotator to pull out HPO terms
    @return - a JSON response of the HPO terms that were identified
    '''
    # annotate with HPO terms using the NCBO annotator
    indications_text = str(request.args.get('indications'))
    payload = {'apikey': '8b5b7825-538d-40e0-9e9e-5ab9274a9aeb',
               'text': indications_text,
               'ontologies': 'HP',
               'whole_word_only': 'false'}
    url = 'http://data.bioontology.org/annotator?' + urlencode(payload, quote_via=quote_plus)
    resp = requests.get(url)

    # get definitions of HPO terms
    hpoTerms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    hpo_terms_defs = {}
    for annot in resp.json():
        hpoid = annot['annotatedClass']['@id'].split('/')[-1].replace('_', ':')
        for hpo, defText in hpoTerms:
            if hpoid.lower() in hpo.lower():
                hpo_terms_defs[hpoid] = defText
                break

    resp.close()

    return jsonify(hpo_terms_defs)

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