#!/usr/bin/env python3

import sys
import argparse
import pyfastaq
import pymummer
import subprocess
import os

parser = argparse.ArgumentParser(
    description = '''Compares FASTA files with blast or nucmer, writes input files for ACT.
                   Then start ACT. Files from top to bottom in ACT are same as order
                   listed on command line when this script is run.''',
    usage = '%(prog)s [options] <blast|nucmer|promer> <outdir> <file1.fa> <file2.fa>  [<file3.fa ...]')
parser.add_argument('--blast_ops', help='blastall options [%(default)s]', default='-p blastn -m 8 -F F -e 0.01 -b 10000 -v 10000')
parser.add_argument('--nucmer_ops', help='nucmer or promer options [promer:--maxmatch. nucmer: --maxmatch --nosimplify]')
parser.add_argument('--no_delta_filter', action='store_true')
parser.add_argument('--no_act', action='store_true', help='Do not start act, just make comparison files etc')
parser.add_argument('--delta_ops', help='delta-filter options [%(default)s]', default='-m')
parser.add_argument('aln_tool', help='blast, nucmer or promer')
parser.add_argument('outdir', help='Output directory (must not already exist)')
parser.add_argument('fa_list', help='List of fasta files', nargs=argparse.REMAINDER)
options = parser.parse_args()
assert len(options.fa_list) > 1


def index_to_union(ops, i):
    return os.path.join(ops.outdir, 'infile.' + str(i) + '.union.fa')


def compare_with_blast(qry, ref, ops, outfile):
    subprocess.check_output('formatdb -l ' + os.path.join(ops.outdir, '.formatdb.log') + ' -p F -i ' + ref, shell=True)
    cmd = ' '.join([
        'blastall', ops.blast_ops,
        '-d', ref,
        '-i', qry,
        '-o', outfile
    ])
    subprocess.check_output(cmd, shell=True)


def compare_with_nucmer(qry, ref, ops, outfile):
    nucmer_out = outfile + '.nucmer.out'
    delta_file = nucmer_out + '.delta'
    filtered_file = delta_file + '.filter'
    coords_file = filtered_file + '.coords'

    if ops.nucmer_ops is None:
        if ops.aln_tool == 'promer':
            ops.nucmer_ops = '--maxmatch'
        else:
            ops.nucmer_ops = '--maxmatch --nosimplify'

    cmd = ' '.join([
        ops.aln_tool,
        ops.nucmer_ops,
        '-p', nucmer_out,
        ref,
        qry,
    ])
    print('cmd:', cmd)
    pyfastaq.utils.syscall(cmd)

    if ops.no_delta_filter:
        cmd = 'cp ' + delta_file + ' ' + filtered_file
    else:
        cmd = ' '.join([
            'delta-filter',
            ops.delta_ops,
            delta_file,
            '>', filtered_file,
        ])

    print('cmd:', cmd)
    pyfastaq.utils.syscall(cmd)

    cmd = ' '.join([
        'show-coords -dTlroH',
        filtered_file,
        '>', coords_file
    ])
    print('cmd:', cmd)
    pyfastaq.utils.syscall(cmd)

    pyfastaq.utils.syscall('samtools faidx ' + qry)
    pyfastaq.utils.syscall('samtools faidx ' + ref)
    pymummer.coords_file.convert_to_msp_crunch(coords_file, outfile, qry + '.fai', ref + '.fai')


# check files exist
for i in range(len(options.fa_list)):
    if not os.path.exists(options.fa_list[i]):
        print('File not found:', options.fa_list[i], file=sys.stderr)
        sys.exit(1)

    options.fa_list[i] = os.path.abspath(options.fa_list[i])

try:
    os.mkdir(options.outdir)
except:
    print('Error making output directory', options.outdir)
    sys.exit(1)


# make union files
for i in range(len(options.fa_list)):
    seq = pyfastaq.sequences.Fasta('union', '')
    reader = pyfastaq.sequences.file_reader(options.fa_list[i])
    new_seq = []
    for s in reader:
        new_seq.append(s.seq)
    f = pyfastaq.utils.open_file_write(index_to_union(options, i))
    seq.seq = ''.join(new_seq)
    print(seq, file=f)
    pyfastaq.utils.close(f)


act_command = 'act ' + options.fa_list[0]

# run alignments
for i in range(len(options.fa_list)-1):
    qry = index_to_union(options, i+1)
    ref = index_to_union(options, i)
    outfile = 'compare.' + str(i) + '.vs.' + str(i+1)
    outfile_abs = os.path.join(options.outdir, outfile)

    if options.aln_tool == 'blast':
        compare_with_blast(qry, ref, options, outfile_abs)
    elif options.aln_tool in ['nucmer', 'promer']:
        compare_with_nucmer(qry, ref, options, outfile_abs)
    else:
        sys.exit('Unknown alignment tool:' + options.aln_tool)

    act_command += ' ' + outfile + ' ' + options.fa_list[i+1]

# delete temporary union files
for i in range(len(options.fa_list)):
    filename = index_to_union(options, i)
    os.unlink(filename)
    try:
        os.unlink(filename + '.fai')
    except:
        pass


# write ACT script
try:
    os.chdir(options.outdir)
except:
    print('Error chdir', options.outdir)
    sys.exit(1)

act_script = 'start_act.sh'
with open(act_script, 'w') as f:
    print('#!/usr/bin/env bash', file=f)
    print('set -e', file=f)
    print('dir=$(dirname $0)', file=f)
    print('cd $dir', file=f)
    print(act_command, file=f)
os.chmod(act_script, 0o755)

if not options.no_act:
    subprocess.check_output('./' + act_script, shell=True)

