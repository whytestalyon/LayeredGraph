'''
This will read the datasets we use and consolidate versions into a single JSON we can load from the server
'''
import json

import pronto

def getOboVersions(ontFN):
    ont = pronto.Ontology(ontFN)
    return (ont.meta['format-version'], ont.meta['data-version'])

def getOwlVersion(owlFN):
    ont = pronto.Ontology(owlFN)
    return ont.meta['versionInfo']

if __name__ == '__main__':
    hpoFN = './HPO_data_files/hp.obo'
    doFN = './HPO_data_files/doid.obo'
    orphaFN = './HPO_data_files/ordo_orphanet.owl'
    
    outFN = './HPO_graph_data/metadata.json'
    
    hpo_fv, hpo_dv = getOboVersions(hpoFN)
    do_fv, do_dv = getOboVersions(doFN)
    orpha_fv = getOwlVersion(orphaFN)
    
    jDict = {
        'hpo-format-version': hpo_fv,
        'hpo-data-version': hpo_dv,
        'do-format-version': do_fv,
        'do-data-version': do_dv,
        'orpha-versionInfo': orpha_fv
    }
    
    fp = open(outFN, 'w+')
    json.dump(jDict, fp)
    fp.close()