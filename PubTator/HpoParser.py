import PubTatorParser


def get_hpo_disease2hpoId_map(source_tags):
    # load hpo disease to phenotype information
    hpo_annot_file = open('./HPO_data_files/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt', 'r')
    hpo_disease_dict = dict({})
    available_tags = source_tags.keys()
    for line in hpo_annot_file:
        cols = line.split("\t")
        recognized_disease_tag = False
        for tag in available_tags:
            if tag in cols[0]:
                recognized_disease_tag = True
                break

        if len(cols) > 3 and cols[3].startswith('HP') and recognized_disease_tag:
            disease_id = PubTatorParser.clean_tag(source_tags, cols[0])
            if disease_id in hpo_disease_dict:
                hpo_disease_dict[disease_id].add(cols[3])
            else:
                hpo_disease_dict[disease_id] = {cols[3]}

    hpo_annot_file.close()

    return hpo_disease_dict
