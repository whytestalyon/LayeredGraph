'''
This will server as a basic wrapper for the LayeredGraph.  It will primarily be an API with minimal actual HTML utility.
'''
import json
from flask import Flask
from flask import g
from flask import jsonify
from flask import render_template
from flask import Response
from flask import request
import pickle
import requests
from urllib.parse import urlencode, quote_plus

# since all creation of Graphs is done from a LayeredGraphAPI subfolder, we have to add it to preserve the pickle data
import sys

sys.path.append("LayeredGraphAPI")
from HPOParser import getTermsAndSyns
from PubTator import HGNCParser
from PhenotypeAPI import PhenotypeCorrelationParser

# define appilication object
app = Flask(__name__)


@app.route('/background')
def background():
    '''
    TODO: Function to return the background array, can be removed later OR folded into the initialization so background isn't constantly recalculated
    '''
    # primary data
    mg = getMultigraphVars()
    hpoWeights = pickle.load(open('./HPO_graph_data/multiHpoWeight_biogrid_pushup.pickle', 'rb'))
    restartProb = 0.1

    # background calculation
    bgNodes = mg.nodes['HPO']
    bgProbs = {('HPO', h): hpoWeights[h] for h in mg.nodes['HPO']}
    bg = mg.calculateBackground(bgProbs, restartProb)

    return jsonify(bg)


@app.route('/')
@app.route('/table', methods=['GET'])
def table():
    '''
    This is the page for the table view, it performs a more complex query so it is given its own page
    '''
    return render_template('textToTable.html')


@app.route('/ppi', methods=['GET'])
def ppi():
    '''
    This will be a deeprank view for the PPI graph.
    '''
    return render_template('protSearch.html')


@app.route('/about', methods=['GET'])
def about():
    '''
    Get the about page
    '''
    return render_template('about.html')


@app.route('/terms', methods=['GET'])
def terms():
    '''
    This function returns a JSON presentation of the terms that match user typed text
    '''
    # pprint(request.args)
    searchTerm = str(request.args.get('term'))
    options = list([])

    for hpo, defText in hpoTerms:
        if searchTerm.lower() in defText.lower():
            options.append({'id': hpo, 'text': (hpo + ' ' + defText)})
        if searchTerm.lower() in hpo.lower():
            options.append({'id': hpo, 'text': (hpo + ' ' + defText)})

    return jsonify({'results': options})


@app.route('/genes', methods=['GET'])
def genes():
    '''
    This function returns a JSON presentation fo the genes that match the supplied search string
    :return:
    '''
    options = list([])
    search_symbol = str(request.args.get('term'))
    for gene_symbol in list(entrez_gene_dict.keys()):
        if search_symbol.lower() in gene_symbol.lower():
            options.append({'id': gene_symbol, 'text': gene_symbol})

    return jsonify({'results': options})


@app.route('/phenotypegene', methods=['POST'])
def phenotypes():
    '''
    Given an input gene and phenotype terms return the publications and phenotypes PyxisMap associates with it
    :return:
    '''
    req_json = request.get_json()
    search_phenotypes_list = req_json['phenotypes']

    phenotypes4gene = PhenotypeCorrelationParser \
        .read_pubmed_info_from_index("./HPO_graph_data/gene2phenotype.json",
                                     gene2phenotype2pub_blockindex[str(entrez_gene_dict[req_json['gene']])],
                                     search_phenotypes_list)

    pmid_dict = dict({})
    for phenotype, pmids in phenotypes4gene.items():
        for pmid in pmids:
            if pmid in pmid_dict:
                pmid_dict[pmid].append(phenotype)
            else:
                pmid_dict[pmid] = [phenotype]

    res = []
    for pmid, phenotypes in pmid_dict.items():
        res.append(
            {'terms': phenotypes, 'pmid': pmid, 'count': len(phenotypes)})

    return jsonify({'data': res})


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

            # do an individual weighting for each term as we go
            termWeights = mg.RWR_rank({('HPO', h): 1.0}, restartProb, rankTypes, bg)
            for i, (w, t, l) in enumerate(termWeights):
                if (l not in indivRanks):
                    indivRanks[l] = {}
                indivRanks[l][h] = (i + 1, w)
        else:
            missingTerms.add(h)

    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)

    rankings = []
    for i, (w, t, l) in enumerate(rankedGenes):
        temp = {'weight': w, 'nodeType': t, 'label': l, 'rank': i + 1}
        temp.update(indivRanks[l])
        rankings.append(temp)

    ret = {'rankings': rankings,
           'usedTerms': list(usedTerms),
           'missingTerms': list(missingTerms)}
    return jsonify(ret)


@app.route('/protrank', methods=['POST'])
def protrank():
    '''
    This is intended to be the workhorse function.  POST needs to include a "term" list that should be HPO terms matching values from the layered graph.
    '''
    mydata = request.get_json()

    # "constant" global data
    mg, restartProb, bg = getProtgraphVars()
    genes = set([str(x) for x in mydata])
    rankTypes = set(['gene'])

    startProbs = {}
    usedTerms = set([])
    missingTerms = set([])
    for gene in genes:
        if gene in mg.nodes['gene']:
            startProbs[('gene', gene)] = 1.0
            usedTerms.add(gene)
        else:
            missingTerms.add(gene)

    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)

    rankings = []
    for w, t, l in rankedGenes:
        rankings.append({'weight': w, 'nodeType': t, 'label': l})

    ret = {'rankings': rankings,
           'usedTerms': list(usedTerms),
           'missingTerms': list(missingTerms)}
    return jsonify(ret)


@app.route('/protdeeprank', methods=['POST'])
def protdeeprank():
    '''
    PPI workhorse, POST needs to include a "term" list that are gene names for the PPI graph.
    '''
    mydata = request.get_json()

    # "constant" global data
    mg, restartProb, bg = getProtgraphVars()
    genes = set([str(x) for x in mydata])
    rankTypes = set(['gene'])

    startProbs = {}
    usedTerms = set([])
    missingTerms = set([])

    indivRanks = {}

    for gene in genes:
        if gene in mg.nodes['gene']:
            startProbs[('gene', gene)] = 1.0
            usedTerms.add(gene)

            # do an individual weighting for each term as we go
            termWeights = mg.RWR_rank({('gene', gene): 1.0}, restartProb, rankTypes, bg)
            for i, (w, t, l) in enumerate(termWeights):
                if (l not in indivRanks):
                    indivRanks[l] = {}
                indivRanks[l][gene] = (i + 1, w)
        else:
            missingTerms.add(gene)

    rankedGenes = mg.RWR_rank(startProbs, restartProb, rankTypes, bg)

    rankings = []
    for i, (w, t, l) in enumerate(rankedGenes):
        temp = {'weight': w, 'nodeType': t, 'label': l, 'rank': i + 1}
        temp.update(indivRanks[l])
        rankings.append(temp)

    graphLimit = 25
    graphNodes = []
    graphEdges = []
    for i, (w, t, l) in enumerate(rankedGenes[:graphLimit]):
        graphNodes.append(l)
        for j, (w2, t2, l2) in enumerate(rankedGenes[:graphLimit]):
            ew = mg.getEdge(t, l, t2, l2, True)
            conf = mg.getEdge(t, l, t2, l2, False)
            if ew > 0:
                graphEdges.append((l, l2, ew, conf))

    ret = {'rankings': rankings,
           'usedTerms': list(usedTerms),
           'missingTerms': list(missingTerms),
           'graphNodes': graphNodes,
           'graphEdges': graphEdges}
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
    try:
        resp = requests.get(url)
    except Exception as e:
        print("Python Exception: " + str(e))
        return jsonify({'annotatorStatus': 'UNKNOWN', 'terms': {}})

    returnCode = resp.status_code

    # get definitions of HPO terms
    hpo_terms_defs = {}
    try:
        for annot in resp.json():
            hpoid = annot['annotatedClass']['@id'].split('/')[-1].replace('_', ':')
            for hpo, defText in hpoTerms:
                if hpoid.lower() in hpo.lower():
                    hpo_terms_defs[hpoid] = defText
                    break
    except Exception as e:
        print("Annotator response: " + str(resp.json()))
        print("Python Exception: " + str(e))

    resp.close()

    return jsonify({'annotatorStatus': returnCode, 'terms': hpo_terms_defs})


def getMultigraphVars():
    '''
    Get any variables that are constant between requests from flask.g.
    Note: I don't think this is doing quite what I want it to because the files are being loaded with each request, but it is relatively fast so moving on for now.
    @return - tuple of variables (multigraph, restart prob, hpo weights, background weights)
    '''
    if not hasattr(g, 'mg'):
        # load or generate the graph
        pickleGraphFN = './HPO_graph_data/multigraph.pickle'
        print('Loading from "' + pickleGraphFN + '"')

        # these are variable that need to be set for multiple uses
        g.mg = pickle.load(open(pickleGraphFN, 'rb'))

    if not hasattr(g, 'restartProb'):
        g.restartProb = 0.1

    if not hasattr(g, 'hpoWeights'):
        g.hpoWeights = pickle.load(open('./HPO_graph_data/multiHpoWeight_biogrid_pushup.pickle', 'rb'))

    if not hasattr(g, 'bg'):
        bgProbs = {('HPO', h): g.hpoWeights[h] for h in g.mg.nodes['HPO']}
        g.bg = g.mg.calculateBackground(bgProbs, g.restartProb)

    return (g.mg, g.restartProb, g.hpoWeights, g.bg)


def getProtgraphVars():
    '''
    Get any variables related to the protein graph from flask.g.
    @return - tuple of variables (multigraph, restart prob, background weights)
    '''
    if not hasattr(g, 'protgraph'):
        # load the graph
        pickleGraphFN = './HPO_graph_data/protgraph.pickle'
        print('Loading from "' + pickleGraphFN + '"')
        g.protgraph = pickle.load(open(pickleGraphFN, 'rb'))

    if not hasattr(g, 'protRestartProb'):
        g.protRestartProb = 0.1

    if not hasattr(g, 'protBG'):
        bgProbs = {('gene', v): 1.0 for v in g.protgraph.nodes['gene']}
        g.protBG = g.protgraph.calculateBackground(bgProbs, g.protRestartProb)

    return (g.protgraph, g.protRestartProb, g.protBG)


# run the app.
if __name__ == "__main__":
    # Setting debug to True enables debug output. This line should be
    # removed before deploying a production app.

    # initializeLayeredGraph()
    # define global data caches
    print('Caching HPO...')
    global hpoTerms
    hpoTerms = getTermsAndSyns('./HPO_graph_data/hp.obo')
    print('Caching HGNC gene data...')
    global entrez_gene_dict
    entrez_gene_dict = HGNCParser.load_genes()
    print('Caching PubTator data...')
    global gene2phenotype2pub_blockindex
    infile = open('./HPO_graph_data/gene2phenotypeIndex.json')
    gene2phenotype2pub_blockindex = json.load(infile)
    infile.close()

    app.run(debug=True, host='0.0.0.0')
