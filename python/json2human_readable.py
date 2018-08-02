#!/usr/bin/env python3

import argparse
import json

parser = argparse.ArgumentParser(
    description = 'Prints human-readable json file, by adding line breaks and indents',
    usage='%(prog)s <in.json>')

parser.add_argument('json_file', help='Name of input json file')
options = parser.parse_args()

with open(options.json_file) as f:
    json_data = json.load(f)
    print(json.dumps(json_data, sort_keys=True, indent=4))

