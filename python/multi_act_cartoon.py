#!/usr/bin/env python3

import os
import sys
import argparse
import pyfastaq
import pymummer


# coords = list of tuples [(x1, y1), (x2, y2) ...]
def svg_polygon(coords, fill_colour, border_colour, border_width = 1, opacity=-1):
    return_string = '<polygon points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
        + 'fill="' + fill_colour + '" '

    if opacity != -1:
        return_string += 'fill-opacity="' + str(opacity) + '" '

    return_string += 'stroke="' + border_colour + '" ' \
                     + 'stroke-width="' + str(border_width) + '" ' \
                     + '/>'
    return return_string


class Appearance:
    def __init__(self, opts):
        self.contig_fill_colour = opts.seq_col
        self.contig_border_colour = opts.seq_col
        self.contig_border_width = 0
        self.contig_x_space = opts.contig_x_space
        self.match_fwd_fill_colour = opts.fwd_match_col
        self.match_fwd_border_colour = opts.fwd_match_col
        self.match_fwd_border_width = 0
        self.match_rev_fill_colour = opts.rev_match_col
        self.match_rev_border_colour = opts.rev_match_col
        self.match_rev_border_width = 0
        self.match_opacity = opts.match_opacity
        self.match_min_length_bases = 10000
        self.match_min_length_ratio = 0.2


class Assembly:
    def __init__(self, fasta_file):
        self.names = []
        self.lengths = {}
        seq_reader = pyfastaq.sequences.file_reader(fasta_file)
        for seq in seq_reader:
            self.lengths[seq.id] = len(seq)
            self.names.append(seq.id)


    def x_width(self, contig_x_space):
        return sum(self.lengths.values()) + contig_x_space * (len(self.lengths) - 1)


    def contigs_svg(self, scale_factor, y_top, y_bottom, appearance):
        x_in_bases = 0
        lines = []
        self.contig_coords = {}

        for contig in self.names:
            x_end = x_in_bases + self.lengths[contig]
            x_left = x_in_bases * scale_factor
            x_right = x_end * scale_factor
            self.contig_coords[contig] = (x_left, x_right)
            coords = [(x_left, y_top), (x_right, y_top), (x_right, y_bottom), (x_left, y_bottom)]
            lines.append(
                svg_polygon(
                    coords,
                    appearance.contig_fill_colour,
                    appearance.contig_border_colour,
                    border_width=appearance.contig_border_width
                )
            )
            x_in_bases = x_end + appearance.contig_x_space

        return '\n'.join(lines)


class Assemblies:
    def __init__(self, fasta_files, appearance, outprefix, nucmer_min_id=98, nucmer_min_length=250, simplify=True):
        self.fasta_files = fasta_files
        self.appearance = appearance
        self.outprefix = outprefix
        ok = True

        for filename in fasta_files:
            if not os.path.exists(filename):
                print('Could not find file:', filename, file=sys.stderr)
                ok = False

        if not ok:
            sys.exit(1)

        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_length = nucmer_min_length
        self.simplify = simplify
        self.assemblies = {filename: Assembly(filename) for filename in fasta_files}


    @staticmethod
    def _get_x_max(assemblies, contig_x_space):
        return max([a.x_width(contig_x_space) for a in assemblies.values()])


    def _make_nucmer_files(self, outprefix):
        self.nucmer_matches = []
        self.nucmer_files = []

        for i in range(len(self.assemblies) - 1):
            print('Comparing assembly ', i, ' (', self.fasta_files[i], ') against ', i+1, ' (', self.fasta_files[i+1], ')', sep='')
            nucmer_file = '.'.join([outprefix, str(i), str(i+1), 'coords'])
            if os.path.exists(nucmer_file):
                print('Found nucmer coords file', nucmer_file, 'so no need to run nucmer')
            else:
                print('Running nucmer. Coords file will be called:', nucmer_file)
                n = pymummer.nucmer.Runner(
                    self.fasta_files[i+1],
                    self.fasta_files[i],
                    nucmer_file,
                    min_id=self.nucmer_min_id,
                    breaklen=500,
                    maxmatch=True,
                    simplify=True,
                    verbose=False,
                )
                n.run()

            self.nucmer_matches.append([x for x in pymummer.coords_file.reader(nucmer_file)])


    @staticmethod
    def _write_svg_header(filehandle, width, height):
        print(r'''<?xml version="1.0" standalone="no"?>''', file=filehandle)
        print(r'''<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">''', file=filehandle)
        print(r'<svg width="' + str(width) + '" height="' + str(height) + '">', file=filehandle)



    def _write_svg_contigs(self, filehandle):
        y_top = 0
        for filename in self.fasta_files:
            y_bottom = y_top + self.contig_height
            print(self.assemblies[filename].contigs_svg(self.x_scale_factor, y_top, y_bottom, self.appearance), file=filehandle)
            y_top += self.contig_height + 2 * self.y_space + self.match_height


    def _write_svg_matches_between_two_assemblies(self, nucmer_matches, top_assembly, bottom_assembly, y_top):
        # top assembly = the query in nucmer matches
        # bottom assembly = the reference in nucmer matches
        lines = []
        y_bottom = y_top + self.match_height
        for match in nucmer_matches:
            if match.hit_length_qry < min(self.appearance.match_min_length_bases, self.appearance.match_min_length_ratio * min(match.qry_length, match.ref_length)):
                continue

            top_contig_start, top_contig_end = top_assembly.contig_coords[match.qry_name]
            top_start = top_contig_start + (match.qry_start / match.qry_length) * (top_contig_end - top_contig_start)
            top_end = top_contig_start + (match.qry_end / match.qry_length) * (top_contig_end - top_contig_start)

            bottom_contig_start, bottom_contig_end = bottom_assembly.contig_coords[match.ref_name]
            bottom_start = bottom_contig_start + (match.ref_start / match.ref_length) * (bottom_contig_end - bottom_contig_start)
            bottom_end = bottom_contig_start + (match.ref_end / match.ref_length) * (bottom_contig_end - bottom_contig_start)

            coords = [(top_start, y_top), (top_end, y_top), (bottom_end, y_bottom), (bottom_start, y_bottom)]
            if match.on_same_strand():
                lines.append((top_start - top_end, svg_polygon(coords, self.appearance.match_fwd_fill_colour, self.appearance.match_fwd_border_colour, opacity=self.appearance.match_opacity, border_width=self.appearance.match_fwd_border_width)))
            else:
                lines.append((top_start - top_end, svg_polygon(coords, self.appearance.match_rev_fill_colour, self.appearance.match_rev_border_colour, opacity=self.appearance.match_opacity, border_width=self.appearance.match_rev_border_width)))

        lines.sort()
        return '\n'.join([x[1] for x in lines])


    def _write_all_svg_matches(self, filehandle):
        y_top = self.contig_height + self.y_space

        for i in range(len(self.fasta_files) - 1):
            top_assembly = self.assemblies[self.fasta_files[i]]
            bottom_assembly = self.assemblies[self.fasta_files[i+1]]
            print(self._write_svg_matches_between_two_assemblies(self.nucmer_matches[i], top_assembly, bottom_assembly, y_top), file=filehandle)
            y_top += self.contig_height + 2 * self.y_space + self.match_height


    def run(self, outprefix):
        self._make_nucmer_files(outprefix)
        self.contig_height = 1.5
        self.match_height = 30
        self.y_space = 1
        self.svg_height = (len(self.assemblies) - 1) * (self.contig_height + 2 * self.y_space + self.match_height) + self.y_space + self.contig_height
        self.svg_width = 400
        self.total_width_in_bases = self._get_x_max(self.assemblies, self.appearance.contig_x_space)
        self.x_scale_factor = self.svg_width / self.total_width_in_bases

        svg_file = outprefix + '.svg'
        print('Writing SVG file', svg_file)
        svg_fh = pyfastaq.utils.open_file_write(svg_file)
        self._write_svg_header(svg_fh, self.svg_width, self.svg_height)
        self._write_svg_contigs(svg_fh)
        self._write_all_svg_matches(svg_fh)
        print('</svg>', file=svg_fh)
        pyfastaq.utils.close(svg_fh)
        print('Finished writing SVG file', svg_file)



parser = argparse.ArgumentParser(
    description = 'Makes cartoon ACT-style figure comparing at least two FASTA files. Files are shown from top to bottom in the same order as listed on the command line when calling this script.',
    usage = '%(prog)s [options] <outprefix> <file1.fa> <file2.fa> [more fasta files ...]')
parser.add_argument('--contig_x_space', type=int, help='Space between each contig, in bases [%(default)s]', default=20000, metavar='INT')
parser.add_argument('--seq_col', help='Colour of sequences [%(default)s]', default='black', metavar='STRING')
parser.add_argument('--fwd_match_col', help='Colour of match on same strands [%(default)s]', default='lightseagreen', metavar='STRING')
parser.add_argument('--rev_match_col', help='Colour of match on opposite strands [%(default)s]', default='peru', metavar='STRING')
parser.add_argument('--match_opacity', type=int, help='Opactiy of matches between 0 and 1. Higher means less transparent [%(default)s]', default=0.8, metavar='FLOAT in [0,1]')
parser.add_argument('--match_min_len_bases', type=int, help='Minimum match length to show. Whether or not a match is shown also depends on the sequence lengths. See --match_min_len_ratio [%(default)s]', default=5000, metavar='INT')
parser.add_argument('--match_min_len_ratio', type=float, help='Minimum match length to show, as proportion of sequence lengths. Using "--match_min_len_bases X --match_min_len_ratio Y" means that a match of length L is shown if L >= min(X, Y * S), where S is max(length of seq1, length of seq2) [%(default)s]', default=0.2, metavar='FLOAT')
parser.add_argument('--nucmer_min_id', type=int, help='Minimum identity when running nucmer [%(default)s]', default=90, metavar='FLOAT')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('fa_list', help='List of at least 2 fasta files', nargs=argparse.REMAINDER)
options = parser.parse_args()


if len(options.fa_list) < 2:
    print('Must have at least two input fasta files! Cannot continue', file=sys.stderr)
    sys.exit(1)

a = Assemblies(
    options.fa_list,
    Appearance(options),
    options.outprefix,
    nucmer_min_id=options.nucmer_min_id,
    simplify=True
)
a.run(options.outprefix)

