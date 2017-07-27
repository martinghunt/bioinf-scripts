#!/usr/bin/env python3

import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description = 'Write BED file of gene coords + a number of upstream bases',
    usage='%(prog)s <gene ids file> <in.genbank> <out.bed>')


parser.add_argument('--add_upstream', type=int, default=0, help='Number of extra bases to add upstream of each gene [%(default)s]')
parser.add_argument('genes_in', help='Name of input file of gene names.')
parser.add_argument('genbank_in', help='Name of input genbank file')
parser.add_argument('outfile', help='Name of output BED file')
options = parser.parse_args()


with open(options.genes_in) as f:
    genes_of_interest = {x.rstrip() for x in f.readlines()}


with open(options.outfile, 'w') as f:
    for seq_record in SeqIO.parse(options.genbank_in, "genbank"):
        for feature in seq_record.features:
            if feature.type == 'gene':
                gene_name = feature.qualifiers.get('gene', [None])[0]
                if gene_name not in genes_of_interest:
                    continue

                strand = {1:'+', -1:'-'}[feature.location.strand]
                if strand == '+':
                    start = max(0, feature.location.start - options.add_upstream)
                    end = feature.location.end
                else:
                    start = feature.location.start
                    end = min(len(seq_record), feature.location.end + options.add_upstream)

                print(seq_record.id, start, end, gene_name, '.', strand, sep='\t', file=f)

