#!/usr/bin/env python3
import argparse
import json
import logging
import requests


def info_from_gi(gi):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    data = {
        "db": "nucleotide",
        "id": gi,
        "rettype": "gb",
        "retmode": "text",
    }

    try:
        r = requests.get(url, data)
    except:
        raise Exception(f"Error getting info from GI {gi}")

    if r.status_code != requests.codes.ok:
        raise Exception(f"Error requesting data. Error code: {r.status_code}. Url: + {r.url}")


    results = {
        "ACCESSION": None,
        "VERSION": None,
        "DBLINK": {},
        "KEYWORDS": [],
        "SOURCE": None,
        "VERSION": None,
    }


    in_dblink = False
    general_names = {"KEYWORDS", "SOURCE", "VERSION"}

    for line in r.text.split("\n"):
        if line.startswith("DBLINK"):
            in_dblink = True
            line = line.split()[1:]
            results["DBLINK"][line[0].rstrip(":")] = " ".join(line[1:])
            continue
        elif in_dblink and line[0] != " ":
            in_dblink = False
        elif in_dblink:
            line = line.split()
            results["DBLINK"][line[0].rstrip(":")] = " ".join(line[1:])
            continue

        fields = line.split()
        if len(fields) < 1:
            continue

        if fields[0] == "ACCESSION":
            results["ACCESSION"] = fields[1:]
        elif fields[0] in general_names:
            line = line.split(maxsplit=1)
            assert len(line) == 2
            results[fields[0]] = line[1]

    return results


def info_from_gis(gi_list):
    results = {}
    for gi in gi_list:
        logging.warning(f"Getting info for {gi}")
        try:
            results[gi] = info_from_gi(gi)
        except:
            results[gi] = None
    return results


parser = argparse.ArgumentParser(
    description="Gets assembly accessions and other info from genbank GI",
    usage="%(prog)s [options] <--gi foo | --gifile FILENAME>",
)
parser.add_argument("--outfmt", choices=["tsv","json"], help="Output format. TSV has essential info, JSON has a bit more [%(default)s]", default="json")
parser.add_argument("--gi", help="GI to look up", metavar="ACCESSION")
parser.add_argument("--gifile", help="File of GIs to look up (one per line)", metavar="FILENAME")

options = parser.parse_args()
if options.gi is None and options.gifile is None:
    raise Exception("Must use --gi or --gifile")
if None not in [options.gi, options.gifile]:
    raise Exception("Do not use both of --gi and --gifile")


if options.gi is not None:
    gi_list = [options.gi]
else:
    with open(options.gifile) as f:
        gi_list = [x.rstrip() for x in f]

results = info_from_gis(gi_list)

if options.outfmt == "json":
    print(json.dumps(results, indent=2))
else:
    print("GI", "ACCESSION", "VERSION", "Assembly", "SOURCE", "KEYWORDS", sep="\t")
    for gi, d in results.items():
        if d["ACCESSION"] is None:
            accession = "None"
        else:
            accession = " ".join(d["ACCESSION"])
        print(
            gi,
            accession,
            d["VERSION"],
            d["DBLINK"].get("Assembly", None),
            d["SOURCE"],
            d["KEYWORDS"],
            sep="\t"
        )
