#!/usr/bin/env python3

import argparse
import json


parser = argparse.ArgumentParser(
    description = 'Gt subset of data from json file, by giving list of keys (e.g. a b c would look for [a][b][c] in json file',
    usage='%(prog)s <in.json> <key1> [key2, key3, ...]')

parser.add_argument('json_in', help='JSON input file')
parser.add_argument('keys', nargs='*', help='List of keys to look for')
options = parser.parse_args()


with open(options.json_in) as f:
    json_data = json.load(f)

d = json_data
for k in options.keys:
    if k in d:
        d = d[k]
    else:
        raise Exception(f'Key not found: {k}')

print(json.dumps(d, sort_keys=True, indent=4))

