import gzip


def get_medgen_disease2hpo():
    medgen_zip_filename = './HPO_data_files/MedGen_HPO_OMIM_Mapping.txt.gz'

    medgendict = dict({})
    with gzip.open(medgen_zip_filename, 'rt') as f_in:
        for line in f_in:
            cols = line.split("|")
            if len(cols) > 5 and cols[5].startswith('HP'):
                key = 'OMIM:' + cols[1]
                if key in medgendict:
                    medgendict[key].add(cols[5])
                else:
                    medgendict[key] = {cols[5]}

    return medgendict
