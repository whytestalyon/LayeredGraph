import glob
import json
import requests

if __name__ == '__main__':
    print('Loading cases...')
    filez = glob.glob('/Users/bwilk/Documents/analysis/PM/UDN/terms/*.txt')

    terms = dict()
    for fp in filez:
        sample_id = fp.split('/')[8].split('_')[0]
        term_list = list()
        terms[sample_id] = term_list
        with open(fp) as file:
            for term in file:
                term_list.append(term.replace('\n', ''))

    print('Getting ranks...')
    url = 'http://pyxis.hudsonalpha.org/rank'
    headers = {'Content-Type': 'application/json'}
    for sid, tl in terms.items():
        if len(tl) == 0:
            print('Skipping ' + sid + ', no terms...')
            continue

        print('Ranking ' + sid + '...')
        try:
            resp = requests.post(url=url, data=json.dumps(tl), headers=headers)
        except Exception as e:
            print("Python Exception: " + str(e))
            exit(1)

        if resp.status_code == 200:
            print('Writing terms to file...')
            with open('/Users/bwilk/Documents/analysis/PM/UDN/ranks/' + sid + '_gene_ranks.txt',
                      'w+') as gene_ranks_file:
                for ranking in resp.json()['rankings']:
                    gene_ranks_file.write(ranking['label'] + '\n')
            resp.close()
        else:
            print('Failed to get terms for: ' + sid + ', resp = ' + str(resp.status_code))
            resp.close()
            exit(2)
