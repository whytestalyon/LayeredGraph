import json


def build_block_index(filepath):
    gene2phenotype_file = open(filepath, "r")

    gene_index = {}
    nextLineByte = gene2phenotype_file.tell()
    line = gene2phenotype_file.readline()
    while line:
        gene2phenotype_json = json.loads(line)
        if gene2phenotype_json["geneId"] not in gene_index:
            gene_index[gene2phenotype_json["geneId"]] = nextLineByte

        nextLineByte = gene2phenotype_file.tell()
        line = gene2phenotype_file.readline()

    gene2phenotype_file.close()

    return gene_index


def read_pubmed_info_from_index(filepath, block_index, accepted_phenotypes):
    phenodict = {}
    with open(filepath, "r") as f:
        f.seek(block_index)
        line = f.readline()
        gene2phenotype_json = json.loads(line)
        gene_of_interest = gene2phenotype_json['geneId']

        f.seek(block_index)
        while line:
            line = f.readline()
            gene2phenotype_json = json.loads(line)
            gene = gene2phenotype_json['geneId']
            phenotype = gene2phenotype_json['hpoId']
            if gene == gene_of_interest:
                if phenotype in accepted_phenotypes:
                    phenodict[gene2phenotype_json["hpoId"]] = gene2phenotype_json["pmids"]
            else:
                break

    return phenodict


if __name__ == "__main__":
    blockindex = build_block_index("/Users/bwilk/workspace/LayeredGraph/HPO_graph_data/gene2phenotype.json")
    outfile = open("/Users/bwilk/workspace/LayeredGraph/HPO_graph_data/gene2phenotypeIndex.json", 'w+')
    json.dump(blockindex, outfile)
    outfile.close()

    outfile = open("/Users/bwilk/workspace/LayeredGraph/HPO_graph_data/gene2phenotypeIndex.json", 'r')
    blockindex = json.load(outfile)
    outfile.close()
    print('block: ' + str(blockindex["79912"]))
    print(str(read_pubmed_info_from_index("/Users/bwilk/workspace/LayeredGraph/HPO_graph_data/gene2phenotype.json",
                                          blockindex["79912"], ["HP:0003198", "HP:0001324"])))
