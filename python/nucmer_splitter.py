#!/usr/bin/env python3

import argparse
import logging
import os
import subprocess


def delta_file_has_a_hit(delta_file):
    with open(delta_file) as f:
        for i, line in enumerate(f):
            if i == 1:
                assert line == "NUCMER\n"
            elif i > 1:
                return True
    return False


parser = argparse.ArgumentParser(
    description="Run nucmer on reference split across multiple fasta files",
    usage="%(prog)s [options] <ref_dir> <query> <outdir>",
)

parser.add_argument(
    "--nuc_opts",
    metavar="STR",
    default="",
    help="delta-filter options string. This is passed straight into nucmer, not sanity checked",
)
parser.add_argument(
    "--delta_opts",
    metavar="STR",
    default="-i 80 -l 100 -m",
    help="delta-filter options string. This is passed straight into delta-filter, not sanity checked [%(default)s]",
)
parser.add_argument(
    "ref_dir",
    help="Directory of reference fasta files. Will use all files in this directory",
)
parser.add_argument("query", help="Query fasta file")
parser.add_argument("outdir", help="Output directory")
options = parser.parse_args()

options.ref_dir = os.path.abspath(options.ref_dir)
options.query = os.path.abspath(options.query)


logging.basicConfig(
    format="[%(asctime)s nucmer_spliiter %(levelname)s] %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
log = logging.getLogger()
log.setLevel(logging.INFO)

if os.path.exists(options.query):
    logging.info(f"Query file: {options.query}")
else:
    raise FileNotFoundError(f"Query file not found {options.query}")

if os.path.exists(options.ref_dir):
    logging.info(f"Reference directory: {options.ref_dir}")
else:
    raise FileNotFoundError(f"Reference directory not found {options.ref_dir}")

ref_files = [os.path.join(options.ref_dir, x) for x in os.listdir(options.ref_dir)]

if len(ref_files) == 0:
    raise Exception(f"No files found in reference directory {options.ref_dir}")

os.mkdir(options.outdir)

if options.query.endswith(".gz"):
    new_query = os.path.join(options.outdir, "qry.fa")
    command = f"gunzip -c {options.query} > {new_query}"
    logging.info(f"Query is gzipped. Extracting: {command}")
    subprocess.check_output(command, shell=True)
    options.query = "qry.fa"
else:
    new_query = None


coords_file = "nucmer.coords"
coords_cols = [
    "[S1]",
    "[E1]",
    "[S2]",
    "[E2]",
    "[LEN 1]",
    "[LEN 2]",
    "[% IDY]",
    "[LEN R]",
    "[LEN Q]",
    "[FRM]",
    "[TAGS]",
    "[NAME R]",
    "[NAME Q]",
    "[EXTRA]",
]
with open(os.path.join(options.outdir, coords_file), "w") as f:
    print(*coords_cols, sep="\t", file=f)

for i, ref in enumerate(ref_files):
    logging.info(f"Processing ref file {i+1}/{len(ref_files)}")
    delta_file = f"{i+1}.delta"
    command = f"nucmer {options.nuc_opts} --delta {delta_file} {ref} {options.query}"
    logging.info(f"  running nucmer: {command}")
    subprocess.check_output(command, shell=True, cwd=options.outdir)

    if not delta_file_has_a_hit(os.path.join(options.outdir, delta_file)):
        logging.info(f"  no hits in {delta_file}. Moving on to next ref file")
        continue

    command = f"delta-filter {options.delta_opts} {delta_file} > {delta_file}.filter"
    logging.info(f"  running delta-filter: {command}")
    subprocess.check_output(command, shell=True, cwd=options.outdir)

    if not delta_file_has_a_hit(os.path.join(options.outdir, f"{delta_file}.filter")):
        logging.info(f"  no hits in {delta_file}.filter. Moving on to next ref file")
        continue

    command = f"show-coords -HdTlro {delta_file}.filter >> {coords_file}"
    logging.info(f"  running show-coords: {command}")
    subprocess.check_output(command, shell=True, cwd=options.outdir)

if new_query is not None:
    logging.info(f"Deleting twmporary unzipped query file {new_query}")
    os.unlink(new_query)


logging.info(f"Finished \N{thumbs up sign}")
