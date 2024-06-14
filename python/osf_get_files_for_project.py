#!/usr/bin/env python3

import argparse
import logging
import requests
import json


def get_one_page_of_data(project_id, path=None, page=1):
    url = f"https://api.osf.io/v2/nodes/{project_id}/files/osfstorage"
    if path is not None:
        url += f"/{path}/"
    opts = {"page[size]": 100, "page": page}
    logging.info(f"Getting url: {url} {opts}")
    response = requests.get(url, opts)
    if response.status_code == 200:
        logging.info(f"Got {response.url}")
        return response.json()
    else:
        raise Exception(f"Error getting url: {response.url}")


def get_all_data(project_id, path=None):
    page = 1
    results = []
    while True:
        new_results = get_one_page_of_data(project_id, path=path, page=page)
        try:
            next_link = new_results["links"]["next"]
        except:
            raise Exception("No 'links/next' entry found in result. Cannot continue")

        results.extend(new_results["data"])
        logging.info(f"New results added. Total results: {len(results)}")

        if next_link is None:
            break
        page += 1

    return results


def get_files_from_folders_iter(project_id, results):
    to_replace = {}
    for i, d in enumerate(results):
        if d["attributes"]["kind"] == "file":
            continue
        elif d["attributes"]["kind"] == "folder":
            path = d["attributes"]["path"].replace("/", "")
            to_replace[i] = get_all_data(project_id, path=path)
        else:
            raise NotImplementedError(
                f"Unknown attribute/kind: {d['attributes']['kind']}\n{d}"
            )

    for i in reversed(sorted(list(to_replace.keys()))):
        results[i : i + 1] = to_replace[i]

    return len(to_replace) > 0


def get_files_from_folders(project_id, results):
    while get_files_from_folders_iter(project_id, results):
        pass


def print_tsv(results):
    lines_out = []
    for d in results:
        lines_out.append(
            (
                d["attributes"]["materialized_path"].lstrip("/"),
                d["links"]["download"],
                d["attributes"]["extra"]["hashes"]["md5"],
                round(float(d["attributes"]["size"]) / (1024 * 1024), 1),
            )
        )
    lines_out.sort()
    print("filename", "url", "md5", "size(MB)", sep="\t")
    for line in lines_out:
        print(*line, sep="\t")


parser = argparse.ArgumentParser(
    description="Get files and metadata for an OSF public project",
    usage="%(prog)s [options] <project_id>",
)
parser.add_argument("project_id", help="Project ID")
parser.add_argument(
    "--format",
    choices=["json", "tsv"],
    default="tsv",
    help="Output format, choose from json (gives all data) or tsv (just outputs filename,url,md5,size) [%(default)s]",
)
options = parser.parse_args()

log = logging.getLogger()
log.setLevel(logging.INFO)

results = get_all_data(options.project_id)
logging.info("Getting folders...")
get_files_from_folders(options.project_id, results)

if options.format == "json":
    print(json.dumps(results, indent=2))
elif options.format == "tsv":
    print_tsv(results)
else:
    raise NotImplementedError(f"Uknown output format {options.format}")
