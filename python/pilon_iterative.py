#!/usr/bin/env python3

import argparse
import errno
import logging
import os
import subprocess


def touch_file(filename):
    with open(filename, 'a'):
        pass


def get_iteration_files(i, options):
    files = {
        'done_file': f'iteration.{i}.done',
        'pilon_dir': f'iteration.{i}.pilon',
        'mapping_prefix': f'iteration.{i}.map',
    }

    if i == 1:
        files['ref_fasta'] = options.assembly_fasta
    else:
        files['ref_fasta'] = os.path.join(f'iteration.{i-1}.pilon', 'pilon.fasta')

    files['corrected_fasta'] = os.path.join(files['pilon_dir'], 'pilon.fasta')
    files['changes_file'] = os.path.join(files['pilon_dir'], 'pilon.changes')
    return files


def log_and_run_command(command):
    logging.info(f'Run: {command}')
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    if completed_process.returncode == 0:
        logging.info(f'Run ok: {command}')
    else:
        logging.error(f'Error running: {command}')
        logging.error(f'Return code: {completed_process.returncode}')
        logging.error(f'Output from stdout:\n{completed_process.stdout}')
        logging.error(f'Output from stderr:\n{completed_process.stderr}')
        raise Exception('Error in system call. Cannot continue')


def make_pilon_bam(reads1, reads2, ref_fasta, outprefix, threads=1):
    sam_file = f'{outprefix}.sam'
    sorted_bam = f'{outprefix}.sorted.bam'
    done_file = f'{outprefix}.done'
    if os.path.exists(done_file):
        logging.info(f'Found done file {done_file}')
        assert os.path.exists(sorted_bam)
        return sorted_bam

    log_and_run_command(f'bwa mem -t {threads} -x intractg {ref_fasta} {reads1} {reads2} > {sam_file}')
    log_and_run_command(f'samtools sort --threads {threads} --reference {ref_fasta} -o {sorted_bam} {sam_file}')
    os.unlink(sam_file)
    log_and_run_command(f'samtools index {sorted_bam}')
    touch_file(done_file)
    return sorted_bam


def run_pilon(bam, ref_fasta, outdir, java_jar, java_xmx_opt, threads=1):
    fasta = os.path.join(outdir, 'pilon.fasta')
    done_file = f'{outdir}.done'
    if os.path.exists(done_file):
        logging.info('Found done file {done_file}')
        assert os.path.exists(f'{fasta}.fai')

    log_and_run_command(f'java -Xmx{java_xmx_opt} -jar {java_jar} --outdir {outdir} --genome {ref_fasta} --frags {bam} --changes --threads {threads} --minmq 10 --minqual 10')
    log_and_run_command(f'bwa index {fasta}')
    log_and_run_command(f'samtools faidx {fasta}')
    touch_file(done_file)


def number_of_pilon_changes(changes_file):
    changes = 0
    with open(changes_file) as f:
        for line in f:
            changes += 1
    return changes


def check_file_exists(to_check, filename_for_log):
    if os.path.exists(to_check):
        logging.info(f'Found {filename_for_log}: {to_check}')
    else:
        logging.info(f'Not found: {filename_for_log}: {to_check}')
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), to_check)


parser = argparse.ArgumentParser(
    description = 'Iteratively run pilon to correct assembly',
    usage='%(prog)s <assembly_fasta> <reads1> <reads2> <outdir>')

parser.add_argument('--pilon_java_xmx', default='8G', help='Option used with -Xmx when running javar pilon.jar [%(default)s]', metavar='STRING')
parser.add_argument('--threads', type=int, default=1, help='Number of threads to use with BWA mem and pilon [%(default)s]', metavar='INT')
parser.add_argument('--max_iterations', type=int, default=10, help='Max number of iterations [%(default)s]', metavar='INT')
parser.add_argument('pilon_jar', help='Name of pilon JAR file')
parser.add_argument('assembly_fasta', help='Name of input assembly fasta file')
parser.add_argument('reads1', help='Name of input reads file 1')
parser.add_argument('reads2', help='Name of input reads file 2')
parser.add_argument('outdir', help='Name of output directory (must not already exist')
options = parser.parse_args()

options.assembly_fasta = os.path.abspath(options.assembly_fasta)
options.reads1 = os.path.abspath(options.reads1)
options.reads2 = os.path.abspath(options.reads2)
options.outdir = os.path.abspath(options.outdir)

if not os.path.exists(options.outdir):
    os.mkdir(options.outdir)

os.chdir(options.outdir)

log = logging.getLogger()
log.setLevel(logging.INFO)
fh = logging.FileHandler('log.txt', mode='a')
log = logging.getLogger()
formatter = logging.Formatter('[pilon_iterative %(asctime)s %(levelname)s] %(message)s', datefmt='%d-%m-%Y %H:%M:%S')
fh.setFormatter(formatter)
log.addHandler(fh)

check_file_exists(options.assembly_fasta, 'assembly_fasta')
check_file_exists(options.reads1, 'reads1')
check_file_exists(options.reads2, 'reads2')


for i in range(1, options.max_iterations + 1, 1):
    logging.info(f'Start iteration {i}')
    files = get_iteration_files(i, options)

    if os.path.exists(files['done_file']):
        logging.info(f'Found iteration done file {files["done_file"]}')
    else:
        bam = make_pilon_bam(options.reads1, options.reads2, files['ref_fasta'], files['mapping_prefix'], threads=options.threads)
        run_pilon(bam, files['ref_fasta'], files['pilon_dir'], options.pilon_jar, options.pilon_java_xmx, threads=options.threads)
        os.unlink(bam)
        os.unlink(f'{bam}.bai')

    number_of_changes = number_of_pilon_changes(files['changes_file'])
    logging.info(f'Number of changes at iteration {i}: {number_of_changes}')
    touch_file(files['done_file'])
    logging.info(f'End iteration {i}')

    if number_of_changes == 0:
        logging.info(f'No changes made in iteration {i}. Stopping')
        break

logging.info('Finished')

