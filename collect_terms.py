import csv
import requests
from urllib.parse import urlencode, quote_plus

if __name__ == '__main__':
    print('Loading cases...')
    with open('/Users/brandon/Documents/PM/UDN/clinical_indications.txt') as summaries_file:
        csv_reader = csv.reader(summaries_file, delimiter='\t')
        summaries_dict = dict(csv_reader)

    print('Querying terms...')
    for case_id, summary in summaries_dict.items():
        print('Encoding ' + case_id + ' summary: ' + summary)
        payload = {'indications': summary}
        url = 'http://pyxis.hudsonalpha.org/text/annotate?' + urlencode(payload, quote_via=quote_plus)
        try:
            print('Querying ' + case_id + '...')
            resp = requests.get(url)
        except Exception as e:
            print("Python Exception: " + str(e))
            exit(1)

        if resp.status_code == 200:
            print('Writing terms to file...')
            with open('/Users/brandon/Documents/PM/UDN/terms/' + case_id + '_phenotypes.txt', 'w+') as terms_file:
                for hpo, desc in resp.json()['terms'].items():
                    terms_file.write(hpo + '\n')
            resp.close()
        else:
            print('Failed to get terms for: ' + case_id + ', resp = ' + str(resp.status_code))
            resp.close()
            exit(2)
