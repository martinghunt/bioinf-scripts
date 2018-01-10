#!/usr/bin/env python3

import argparse
import sys

import pysam
import pyfastaq

class Error (Exception): pass

def decode2(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def load_reads(infile):
    reads = {}

    sam_reader = pysam.Samfile(infile, "rb")
    for sam in sam_reader.fetch(until_eof=True):
        if sam.qname not in reads:
            reads[sam.qname] = [{}, {}]

        if sam.is_read1:
            d = reads[sam.qname][0]
        elif sam.is_read2:
            d = reads[sam.qname][1]
        else:
            print('Read', sam.name,  'must be first or second of pair according to flag. Cannot continue', file=sys.stderr)
            sys.exit(1)

        seq = pyfastaq.sequences.Fastq('x', decode2(sam.seq), decode2(sam.qual))
        if sam.is_reverse:
            seq.revcomp()

        key = (seq.seq, seq.qual)
        d[key] = d.get(key, 0) + 1

    return reads


def run(options):
    reads = load_reads(options.bam)
    f1 = pyfastaq.utils.open_file_write(options.fq_out1)
    f2 = pyfastaq.utils.open_file_write(options.fq_out2)
    counts = {}

    for name_prefix, (reads1, reads2) in reads.items():
        if len(reads1) != 1 or len(reads2) != 1:
            print('Wrong number of reads for name', name_prefix, file=sys.stderr)
            print(reads1, reads2, file=sys.stderr)
            pyfastaq.utils.close(f1)
            pyfastaq.utils.close(f2)
            raise Error('Wrong number of reads. Cannot continue')

        read1 = list(reads1.keys())[0]
        read2 = list(reads2.keys())[0]
        counts[reads1[read1]] = counts.get(reads1[read1], 0) + 1
        counts[reads2[read2]] = counts.get(reads2[read2], 0) + 1

        print('@' + name_prefix + '/1', read1[0], '+', read1[1], sep='\n', file=f1)
        print('@' + name_prefix + '/2', read2[0], '+', read2[1], sep='\n', file=f2)

    pyfastaq.utils.close(f1)
    pyfastaq.utils.close(f2)
    with open(options.count_out, 'w') as f:
        for k, v in sorted(counts.items()):
            print(k, v, sep='\t', file=f)


parser = argparse.ArgumentParser(
    description = 'Converts BAM to apired FASTQ. But removes duplicated reads. This is to fix incorrect BAM files, where the same read appears in the BAM more than once. This is not duplicate in the usual sense of samtools rmdup. Assumes no sort order in the BAM and loads all reads into memory, so can be memory hungry!',
    usage='%(prog)s <in.bam> <out_1.fq> <out_2.fq> <out.count.tsv',
)
parser.add_argument('bam', help='Input BAM filename')
parser.add_argument('fq_out1', help='Output fwd reads file')
parser.add_argument('fq_out2', help='Output rev reads file')
parser.add_argument('count_out', help='Output counts TSV file')

if __name__ == "__main__":
    options = parser.parse_args()
    run(options)
