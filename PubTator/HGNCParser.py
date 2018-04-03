import re


def get_hgnc_genes_ids():
    hgnc_filename = './HPO_data_files/non_alt_loci_set.txt'

    # load hpo disease to phenotype information
    hgnc_file = open(hgnc_filename, 'r')
    accepted_entrez_ids = set({})
    for line in hgnc_file:
        cols = line.split("\t")
        if len(cols) > 18 and re.match(r'\d', cols[18]):
            accepted_entrez_ids.add(cols[18])

    hgnc_file.close()

    return accepted_entrez_ids

def load_genes():
    hgnc_filename = './HPO_graph_data/non_alt_loci_set.txt'
    hgnc_file = open(hgnc_filename, "r")
    genes_dict = {}
    cntr = 0
    for line in hgnc_file:
        cntr += 1
        if cntr == 1:
            continue
        cols = line.split("\t")
        genes_dict[str(cols[1]).upper()] = cols[18]

    hgnc_file.close()
    return genes_dict