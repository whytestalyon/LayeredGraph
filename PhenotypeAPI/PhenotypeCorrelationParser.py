import json


def get_index():
    # since this file is so huge loading it into memory is not an option, so we are gonna get clever and index the file!
    gene2phenotype_file = open("./HPO_graph_data/gene2phenotype.json", "r")

    gene_index = {}
    nextLineByte = gene2phenotype_file.tell()
    line = gene2phenotype_file.readline()
    while line:
        gene2phenotype_json = json.loads(line)
        if gene2phenotype_json["geneId"] in gene_index:
            gene_index[gene2phenotype_json["geneId"]].append(nextLineByte)
        else:
            gene_index[gene2phenotype_json["geneId"]] = [nextLineByte]

        nextLineByte = gene2phenotype_file.tell()
        line = gene2phenotype_file.readline()

    gene2phenotype_file.close()

    return gene_index


def read_json_from_index_list(indexes):
    phenodict = {}
    with open("./HPO_graph_data/gene2phenotype.json", "r") as f:
        for index in indexes:
            f.seek(index)
            line = f.readline()
            gene2phenotype_json = json.loads(line)
            phenodict[gene2phenotype_json["hpoId"]] = gene2phenotype_json["pmids"]

    return phenodict


if __name__ == "__main__":
    print(get_index("/Users/bwilk/Documents/PyxisMap/l2g/data/gene2phenotype.json")["9693"])
