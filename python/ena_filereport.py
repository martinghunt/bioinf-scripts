#!/usr/bin/env python3

import argparse
import requests
import sys

parser = argparse.ArgumentParser(
    description = 'Queries ENA using their filereport REST API. See https://www.ebi.ac.uk/ena/browse/file-reports',
    usage='%(prog)s <accession> <read_run|analysis> <field1[,field2,...]>',
    epilog='For example, get all sample and run accessions for the project PRJEB5162: %(prog)s PRJEB5162 read_run sample_accession,run_accession',
)

allowed_results = ['read_run', 'analysis']
parser.add_argument('accession', help='Accession ID')
parser.add_argument('result', choices=allowed_results, help='Result type. Must be one of: ' + ','.join(allowed_results))
parser.add_argument('fields', help='comma-separated list of fields. eg run_accession,sample_accession')
options = parser.parse_args()

url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?'
data = {
    'accession': options.accession,
    'result': options.result,
    'download': 'txt',
    'fields': options.fields,
}

try:
    r = requests.get(url, data)
except:
    print('Error querying ENA using data:', data, file=sys.stderr)
    print('URL: ', r.url, file=sys.stderr)
    sys.exit(1)

if r.status_code != requests.codes.ok:
    print('Error requesting data. Error code: ', r.status_code, file=sys.stderr)
    print('URL: ', r.url, file=sys.stderr)
    sys.exit(1)

lines = r.text.rstrip().split('\n')
print(*lines, sep='\n')

