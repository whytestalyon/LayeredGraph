import pronto
import zipfile
import gzip
import json
import networkx as nx
import time
import re

import MedGenParser
import HpoParser
import HGNCParser

# since this the Phenotype API is in a different folder we need to add it to the python path
import sys
sys.path.insert(0, '../PhenotypeAPI/')

import PhenotypeCorrelationParser


def put2dict_of_sets(dict, key, item):
    if key in dict:
        dict[key].add(item)
    else:
        dict[key] = {item}


def get_mesh_ids(*args):
    mesh_set = set({})
    for xref_map in args:
        for xref_set in xref_map.values():
            for xref in xref_set:
                if xref.startswith('MESH'):
                    mesh_set.add(xref)

    return mesh_set


def clean_tag(source_tags_dict, id):
    formatted_id = id
    for source_tag in source_tags_dict.keys():
        if source_tag in id:
            formatted_id = id.replace(source_tag, source_tags_dict[source_tag])
            break

    return formatted_id


def get_xrefs_from_ontology(ontology, source_tags_dict):
    print("Getting xrefs for ID for " + str(ontology))
    term_xref_dict = dict({})
    for term in ontology:
        if 'xref' in term.other:
            cleaned_id = clean_tag(source_tags_dict, term.id)
            xrefs = set({})
            for xref in set(term.other['xref']):
                xrefs.add(clean_tag(source_tags_dict, xref))

            term_xref_dict[cleaned_id] = xrefs

    return term_xref_dict


def get_mesh_xref_map(*args):
    mesh_dict = dict({})
    for xref_map in args:
        for xref_key, xref_set in xref_map.items():
            for xref in xref_set:
                if xref.startswith('MESH'):
                    if xref in mesh_dict:
                        mesh_dict[xref].add(xref_key)
                        mesh_dict[xref].update(xref_set)
                    else:
                        mesh_dict[xref] = {xref_key}
                        mesh_dict[xref].update(xref_set)

    return mesh_dict


if __name__ == '__main__':
    start = time.time()
    # load ontologies used for processing relationships between terms
    print('Loading HPO...')
    hpo_file_name = './HPO_data_files/hp.obo'
    hp_ontology = pronto.Ontology(hpo_file_name)

    print('Loading DO...')
    do_file_name = './HPO_data_files/doid.obo'
    disease_ontology = pronto.Ontology(do_file_name)

    print("Loading ORDO...")
    ordo_owl_file = './HPO_data_files/ordo_orphanet.owl'
    ordo_zip_file_name = './HPO_data_files/ordo_orphanet.owl.zip'
    ordo_zip = zipfile.ZipFile(ordo_zip_file_name)
    ordo_zip.extract('ordo_orphanet.owl', './HPO_data_files/')
    ordo_zip.close()

    orphanet_ontology = pronto.Ontology(ordo_owl_file)

    # list out the tags that can be used for cross referencing
    source_tags = {
        "MSH": "MESH",
        "MeSH": "MESH",
        "ORPHA:": "Orphanet_",
        "OMIM": "OMIM",
        "MESH": "MESH",
        "ORPHANET:": "Orphanet_",
        "Orphanet:": "Orphanet_",
        "ORDO:": "Orphanet_",
        "DOID": "DOID",
        "HP:": "HP:"
    }

    accepted_clean_source_tags = {
        "MESH",
        "Orphanet_",
        "OMIM",
        "DOID",
        "HP:"
    }

    # create maps from a given ID in an ontology to it's cross referenced ID in other classification ontologies
    print('Creating cross reference maps...')
    orpha_xref_map = get_xrefs_from_ontology(orphanet_ontology, source_tags)
    hpo_xref_map = get_xrefs_from_ontology(hp_ontology, source_tags)
    do_xref_map = get_xrefs_from_ontology(disease_ontology, source_tags)

    # create mapping of MESH IDs to other IDs
    mesh2other_map = get_mesh_xref_map(orpha_xref_map, hpo_xref_map, do_xref_map)

    # load HPO annotations of disease to phenotype mappings
    print('Loading HPO disease to phenotype mappings...')
    hpo_disease2hpo_map = HpoParser.get_hpo_disease2hpoId_map(source_tags)

    # load MedGen OMIM disease to phenotype mappings
    print('Loading MedGen data...')
    medgen_disease2hpo_map = MedGenParser.get_medgen_disease2hpo()

    # create a set of MESH terms that we can use for mapping
    print('Getting set of accepted mesh terms...')
    accepted_mesh_terms = get_mesh_ids(orpha_xref_map, hpo_xref_map, do_xref_map)
    print('Number of distinct accepted mesh terms: ' + str(len(accepted_mesh_terms)))

    # create graphs per mesh term to map to phenotype term
    print('Building graphs...')
    mesh2disease_phenotype_map = dict({})
    for mesh_id in accepted_mesh_terms:
        # create one graph per MESH ID term
        mesh2disease_phenotype_map[mesh_id] = nx.Graph()
        # link mesh to accepted other terms
        for id in mesh2other_map[mesh_id]:
            for accepted_clean_source_tag in accepted_clean_source_tags:
                if id.startswith(accepted_clean_source_tag):
                    mesh2disease_phenotype_map[mesh_id].add_edge(mesh_id, id)

        # link other IDs to each other
        id_set = set({})
        for id in list(mesh2disease_phenotype_map[mesh_id].nodes):
            if id in do_xref_map:
                id_set = do_xref_map[id]
            elif id in hpo_xref_map:
                id_set = hpo_xref_map[id]
            elif id in orpha_xref_map:
                id_set = orpha_xref_map[id]

            if len(id_set) > 0:
                for link_id in id_set:
                    for accepted_clean_source_tag in accepted_clean_source_tags:
                        if link_id.startswith(accepted_clean_source_tag):
                            mesh2disease_phenotype_map[mesh_id].add_edge(link_id, id)

        # associate phenotypes with diseases in the graph
        id_set = set({})
        for id in list(mesh2disease_phenotype_map[mesh_id].nodes):
            if not id.startswith('HP'):
                id_set.add(id)

        for disease_id in id_set:
            if disease_id in medgen_disease2hpo_map:
                for phenotype in medgen_disease2hpo_map[disease_id]:
                    mesh2disease_phenotype_map[mesh_id].add_edge(disease_id, phenotype)
            if disease_id in hpo_disease2hpo_map:
                for phenotype in hpo_disease2hpo_map[disease_id]:
                    mesh2disease_phenotype_map[mesh_id].add_edge(disease_id, phenotype)

        if mesh_id == 'MESH:D009135' or mesh_id == 'MESH:D018908' or mesh_id == 'MESH:C567499':
            print('Graph for ' + mesh_id)
            print('graph ' + mesh_id.replace(':', '_') + '{')
            for edge in mesh2disease_phenotype_map[mesh_id].edges:
                print('  ' + edge[0].replace(':', '_') + ' -- ' + edge[1].replace(':', '_') + ';')
            print('}')
            print('')

    # load gene information from HGNC
    accepted_entrez_ids = HGNCParser.get_hgnc_genes_ids()
    print("# accepted entrez gene IDs: " + str(len(accepted_entrez_ids)))

    # create map of all accepted genes to all accepted terms

    # load pubtator data
    print('Loading and processing Pubtator data...')
    gene_filename = './HPO_data_files/gene2pubtator.gz'
    disease_filename = './HPO_data_files/disease2pubtator.gz'
    gene2pubmed = {}
    pubmed2disease = {}

    print('Loading gene-pubmed...')
    with gzip.open(gene_filename, 'rt') as f_in:
        for line in f_in:
            cols = line.rstrip().split('\t')

            pmid = cols[0]
            pubmed2disease[pmid] = set([])

            genes = cols[1]
            for gene in re.split(',|;', genes):
                if (gene in accepted_entrez_ids):
                    put2dict_of_sets(gene2pubmed, gene, pmid)

    print('Loading pubmed-disease...')
    with gzip.open(disease_filename, 'rt') as f_in:
        for line in f_in:
            cols = line.rstrip().split('\t')

            pmid = cols[0]
            if (pmid in pubmed2disease):
                disease = cols[1]
                if (disease in accepted_mesh_terms):
                    put2dict_of_sets(pubmed2disease, pmid, disease)

    print('Loaded', len(gene2pubmed), 'gene-pubmed keys.')
    print('Loaded', len(pubmed2disease), 'pubmed-disease keys.')

    print('Writing results to file...')
    outfile = open("./HPO_data_files/gene2phenotype.json", 'w+')

    for gene in sorted(gene2pubmed.keys()):
        hpo2pubmed = {}
        for pubmed in gene2pubmed[gene]:
            for disease in pubmed2disease.get(pubmed, []):
                for phen_id in mesh2disease_phenotype_map[disease].nodes:
                    if str(phen_id).startswith('HP'):
                        put2dict_of_sets(hpo2pubmed, phen_id, pubmed)

        for hpo in sorted(hpo2pubmed.keys()):
            outfile.write(json.dumps({'geneId': gene, 'hpoId': hpo, 'pmids': sorted(hpo2pubmed[hpo])}) + '\n')

    outfile.close()

    print('Building block index for JSON file...')
    block_index = PhenotypeCorrelationParser.build_block_index("./HPO_data_files/gene2phenotype.json")

    print('Writing block index to file...')
    outfile = open("./HPO_data_files/gene2phenotypeIndex.json", 'w+')
    json.dump(block_index, outfile)
    outfile.close()

    print('All PubTator information consumed in ' + str((time.time() - start)) + ' seconds')
