#!/usr/bin/env python3

import argparse
import pymummer
import sys

parser = argparse.ArgumentParser(
    description = 'Convert mummer coords file to mspcrunch format',
    usage = '%(prog)s [options] <infile> <outfile>')
parser.add_argument('--ref_fai', help='Name of reference fai file. Must be used with --qry_fai. If used, makes ACT-compatible crunch file with coords offset correctly' )
parser.add_argument('--qry_fai', help='Name of query fai file. Must be used with --ref_fai. If used, makes ACT-compatible crunch file with coords offset correctly' )
parser.add_argument('infile')
parser.add_argument('outfile')
options = parser.parse_args()


pymummer.coords_file.convert_to_msp_crunch(
    options.infile,
    options.outfile,
    ref_fai=options.ref_fai,
    qry_fai=options.qry_fai
)

