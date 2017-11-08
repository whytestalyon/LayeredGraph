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
    terms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    return render_template('search.html', terms=terms)

@app.route('/')
@app.route('/text', methods=['GET'])
def text():
    return render_template('textToRank.html')

@app.route('/terms', methods=['GET'])
def terms():
    pprint(request.args)
    searchTerm = str(request.args.get('term'))
    hpoTerms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    options = list([])

    for hpo, defText in hpoTerms:
        if searchTerm.lower() in defText.lower():
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

@app.route('/text/annotate', methods=['GET'])
def textannotate():
    # annotate with HPO terms using the NCBO annotator
    indications_text = str(request.args.get('indications'))
    payload = {'apikey': '8b5b7825-538d-40e0-9e9e-5ab9274a9aeb',
               'text': indications_text,
               'ontologies': 'HP',
               'whole_word_only': 'false'}
    url = 'http://data.bioontology.org/annotator?' + urlencode(payload, quote_via=quote_plus)
    resp = requests.get(url)

    # get definitions of HPO terms from PURL for the term
    hpo_terms_defs = {}
    for annot in resp.json():
        purl = annot['annotatedClass']['@id']
        hpoid = annot['annotatedClass']['@id'].split('/')[-1].replace('_', ':')
        presp = requests.get(purl)
        hpo_info = xmltodict.parse(presp.content)

        for rdf_class in hpo_info['rdf:RDF']['Class']:
            if rdf_class['@rdf:about'] == purl:
                hpo_terms_defs[hpoid] = rdf_class['rdfs:label']['#text']
                break

        presp.close()

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