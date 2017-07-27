#!/usr/bin/env python3

import os
import copy
import sys
import argparse
from operator import itemgetter, attrgetter

import pyfastaq
import vcf
from Bio import SeqIO


class Gene:
    def __init__(self, name, strand, start, end, upstream_start, upstream_end):
        self.name = name
        self.strand = strand
        self.start = start
        self.end = end
        self.upstream_start = upstream_start
        self.upstream_end = upstream_end
        assert self.strand in {-1, 1}

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        return '\t'.join([self.name, str(self.strand), str(self.start + 1), str(self.end + 1)])


def load_gene_coords_and_reference(gene_list_file, genbank_file, add_upstream):
    with open(gene_list_file) as f:
        genes_of_interest = {x.rstrip() for x in f.readlines()}

    ref_genes = {}
    ref_seqs = {}

    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        ref_seqs[seq_record.id] = pyfastaq.sequences.Fasta(seq_record.id, str(seq_record.seq))
        gene_list = []

        for feature in seq_record.features:
            if feature.type == 'gene':
                gene_name = feature.qualifiers.get('gene', [None])[0]
                if gene_name not in genes_of_interest:
                    continue

                if feature.location.strand == 1:
                    start = max(0, feature.location.start - add_upstream)
                    end = int(feature.location.end) - 1
                else:
                    start = int(feature.location.start)
                    end = min(len(seq_record) - 1, feature.location.end + add_upstream - 1)

                gene_list.append(Gene(gene_name, feature.location.strand, int(feature.location.start), int(feature.location.end) - 1, start, end))

        gene_list.sort(key=attrgetter('upstream_start'))
        ref_genes[seq_record.id] = gene_list

    return ref_genes, ref_seqs


def load_vcf_file(vcf_file):
    '''Returns tuple: (sample name, vcf_dict), where vcf_dict keys are reference
    sequence names, and each value is a list of vcf records, sorted by chromosome position'''
    vcf_data = {}
    sample_set = set()
    f = open(vcf_file)
    vcf_reader = vcf.Reader(f)

    for record in vcf_reader:
        if record.INFO['SVTYPE'] not in {'INS', 'DEL', 'COMPLEX'}:
            continue

        assert len(record.samples) == 1
        sample_set.add(record.samples[0].sample)

        if record.CHROM not in vcf_data:
            vcf_data[record.CHROM] = []

        vcf_data[record.CHROM].append(copy.copy(record))

    f.close()

    for ref in vcf_data:
        vcf_data[ref].sort(key=lambda x: x.POS)

    if len(sample_set) != 1:
        print('Error! Should only find one sample name! Got: ' + str(sample_set))
    sample = sample_set.pop()
    return sample, vcf_data


def indels_for_one_gene(ref_name, gene, vcf_data, indels):
    if ref_name not in vcf_data:
        return {}
    gene_interval = pyfastaq.intervals.Interval(gene.upstream_start, gene.upstream_end)

    for vcf_index, vcf_record in enumerate(vcf_data[ref_name]):
        variant_start = vcf_record.POS - 1
        variant_end = variant_start + len(vcf_record.REF) - 1
        variant_interval = pyfastaq.intervals.Interval(variant_start, variant_end)

        if gene_interval.intersects(variant_interval):
            assert gene_interval.intersects(variant_interval)
            if (gene.strand == 1 and variant_start == gene.upstream_end): # or (gene.strand == -1 and variant_end == gene.upstream_start):
                continue
            upstream = variant_end < gene.start or variant_start > gene.end
            if (ref_name, vcf_index) not in indels:
                indels[(ref_name, vcf_index)] = []

            indels[(ref_name, vcf_index)].append((upstream, gene))
        elif variant_interval < gene_interval:
            continue
        elif variant_interval > gene_interval:
            break



def get_all_indels(ref_genes, vcf_data):
    # We need to take care of when an indel is in a gene, but is also
    # upstream of another gene. In which case, don't count it twice.
    # Want to store a dict: vcf record -> set of genes. But vcf records are
    # not hashable, so instead have key=(ref name, index) of vcf record, which is the
    # record in the dict vcf_data[ref_name][index]
    indels_dict = {}

    for ref_name in ref_genes:
        for gene in ref_genes[ref_name]:
            indels_for_one_gene(ref_name, gene, vcf_data, indels_dict)

    # In the case where an indel is in one or more gene(s), but also upstream of
    # another gene, discard the upstream matches
    for key in indels_dict:
        upstream_indexes = {i for i in range(len(indels_dict[key])) if indels_dict[key][i][0]}

        if len(upstream_indexes) < len(indels_dict[key]):
            indels_dict[key] = [indels_dict[key][i] for i in range(len(indels_dict[key])) if i not in upstream_indexes]

        #if len(indels_dict[key]) > len(upstream_keys):
        #    for x in upstream_keys:
        #        indels_dict

        #if len(indels_dict[key]) > 1:
        #    indels_dict[key] = {x for x in indels_dict[key] if not x[0]}

        #assert len(indels_dict[key]) == 1 # otherwise it was in two genes?!
        #indels_dict[key] = indels_dict[key].pop()

    return indels_dict


def revcomp(string):
    seq = pyfastaq.sequences.Fasta('x', string)
    seq.revcomp()
    return seq.seq


def genome_coord_to_gene_coord(genome_coord, gene, var_type):
    if gene.strand == 1:
        gene_coord = genome_coord - gene.start
        if var_type == 'INS':
            gene_coord -= 1
    else:
        gene_coord = gene.end - genome_coord
        if var_type == 'INS':
            gene_coord += 1

    if gene_coord >= 0:
        gene_coord += 1


    return gene_coord



def write_output_file(sample, vcf_data, indels, ref_genome, outfile):
    lines_out = []

    for (ref_name, vcf_index), gene_list in indels.items():
        for (gene_is_upstream, gene) in gene_list:
            vcf_record = vcf_data[ref_name][vcf_index]
            variant_position = vcf_record.POS - 1
            ref_seq = ref_genome[vcf_record.CHROM]
            if vcf_record.REF != ref_seq[variant_position:variant_position + len(vcf_record.REF)]:
                print('Mismatch between VCF REF and the actual reference sequence. SKIPPING', sample, vcf_record.CHROM, vcf_record.POS, vcf_record.REF, sep='\t', file=sys.stderr)
                continue
            assert vcf_record.REF == ref_seq[variant_position:variant_position + len(vcf_record.REF)]
            assert vcf_record.INFO['SVTYPE'] in {'INS', 'DEL', 'COMPLEX'}


            start_in_genome = variant_position + 1
            ref_nucleotides = vcf_record.REF[1:]
            alt_nucleotides = str(vcf_record.ALT[0])[1:]
            if vcf_record.INFO['SVTYPE'] == 'INS':
                end_in_genome = start_in_genome + len(vcf_record.REF)
                mutation_var_string = str(vcf_record.ALT[0])[1:]
            else:
                end_in_genome = start_in_genome + len(vcf_record.REF) - 2
                mutation_var_string = vcf_record.REF[1:]

            start_in_gene = genome_coord_to_gene_coord(start_in_genome, gene, vcf_record.INFO['SVTYPE'])
            end_in_gene = genome_coord_to_gene_coord(end_in_genome, gene, vcf_record.INFO['SVTYPE'])

            if gene.strand == -1:
                start_in_gene, end_in_gene = end_in_gene, start_in_gene
                mutation_var_string = revcomp(mutation_var_string)
                ref_nucleotides = revcomp(ref_nucleotides)
                alt_nucleotides = revcomp(alt_nucleotides)

            if end_in_gene > len(gene):
                end_in_gene = str(len(gene)) + '+' + str(end_in_gene - len(gene))

            if vcf_record.INFO['SVTYPE'] == 'COMPLEX':
                mutation = gene.name + '_' + str(start_in_gene) + '_' + str(end_in_gene) + '_' + 'del' + ref_nucleotides + 'ins' + alt_nucleotides
            else:
                mutation = gene.name + '_' + str(start_in_gene) + '_' + str(end_in_gene) + '_' + vcf_record.INFO['SVTYPE'].lower() + '_' + mutation_var_string
            cov1, cov2 = vcf_record.samples[0].data.COV
            cov_string = ','.join([str(x) for x in vcf_record.samples[0].data.COV])
            filter = ','.join(vcf_record.FILTER) if len(vcf_record.FILTER) else 'PASS'
            lines_out.append((sample, vcf_record.CHROM, variant_position + 1, vcf_record.REF, str(vcf_record.ALT[0]), mutation, gene.strand, vcf_record.samples[0].data.GT_CONF, cov1, cov2, vcf_record.samples[0].data.GT, filter))

    lines_out.sort(key=itemgetter(1, 2))

    f = pyfastaq.utils.open_file_write(outfile)
    print('Sample', 'CHROM', 'POS', 'REF', 'ALT', 'variant', 'gene_strand', 'GT_CONF', 'COV1', 'COV2', 'GT', 'FILTER', sep='\t', file=f)
    for line in lines_out:
        print(*line, sep='\t', file=f)
    pyfastaq.utils.close(f)


def run(options):
    ref_genes, ref_genome = load_gene_coords_and_reference(options.genes_in, options.genbank_in, options.add_upstream)
    sample, vcf_data = load_vcf_file(options.vcf_in)
    indels = get_all_indels(ref_genes, vcf_data)
    write_output_file(sample, vcf_data, indels, ref_genome, options.outfile)


parser = argparse.ArgumentParser(
    description = 'Convert cortex vcf to file of indels for genes of interest',
    usage='%(prog)s <gene ids file> <ref.genbank> <in.vcf> <outfile>')

parser.add_argument('--add_upstream', type=int, default=0, help='Number of extra bases to add upstream of each gene [%(default)s]')
parser.add_argument('genes_in', help='Name of input file of gene names. One line per gene')
parser.add_argument('genbank_in', help='Name of input genbank file')
parser.add_argument('vcf_in', help='Name of input vcf file')
parser.add_argument('outfile', help='Name of output file')

if __name__ == "__main__":
    options = parser.parse_args()
    run(options)
