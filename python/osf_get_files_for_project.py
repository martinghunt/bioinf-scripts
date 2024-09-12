#!/usr/bin/env python3

import argparse
import logging
import requests
import json


def get_project_title(node_id):
    url = f"https://api.osf.io/v2/nodes/{node_id}/"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
    else:
        raise Exception(
            f"Failed to fetch project title for node {node_id}: {response.status_code}"
        )

    try:
        return data["data"]["attributes"]["title"]
    except:
        raise Exception(f"No title found for {node_id}\n{data}")


def get_all_children(node_id):
    url = f"https://api.osf.io/v2/nodes/{node_id}/children/"
    all_nodes = {}

    while url:
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()

            current_level_nodes = {
                child["id"]: {
                    "title": child["attributes"]["title"],
                    "parent_id": node_id,
                }
                for child in data["data"]
            }
            all_nodes.update(current_level_nodes)

            # Recursively get children of the current children
            for child_node_id in current_level_nodes:
                all_nodes.update(get_all_children(child_node_id))

            url = data["links"].get("next")
        else:
            raise Exception(
                f"Failed to fetch components for node {node_id}: {response.status_code}"
            )

    return all_nodes


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


def get_all_files_one_project(project_id, path=None):
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
            to_replace[i] = get_all_files_one_project(project_id, path=path)
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


def print_tsv(results, size_in_bytes=False):
    lines_out = []
    for node_id in results:
        if len(results[node_id]["files"]) == 0:
            continue

        assert results[node_id]["title"] is not None
        names = [results[node_id]["title"]]
        parent = results[node_id]["parent_id"]
        while parent is not None:
            names.append(results[parent]["title"])
            assert names[-1] is not None
            parent = results[parent]["parent_id"]
        names.reverse()
        project_names = "/".join(names)

        for d in results[node_id]["files"]:
            if size_in_bytes:
                size = d["attributes"]["size"]
            else:
                size = round(float(d["attributes"]["size"]) / (1024 * 1024), 2)
            lines_out.append(
                (
                    project_names,
                    node_id,
                    d["attributes"]["materialized_path"].lstrip("/"),
                    d["links"]["download"],
                    d["attributes"]["extra"]["hashes"]["md5"],
                    size,
                )
            )
    lines_out.sort()
    size_name = "size(B)" if size_in_bytes else "size(MB)"
    print("project", "project_id", "filename", "url", "md5", size_name, sep="\t")
    for line in lines_out:
        print(*line, sep="\t")


def get_all_data(project_id):
    logging.info(f"Start collecting data for {project_id}")
    data = get_all_children(project_id)
    logging.info(
        f"Found all descendents of {options.project_id}: "
        + ",".join(sorted(list(data.keys())))
    )
    assert project_id not in data
    data[project_id] = {
        "parent_id": None,
        "title": get_project_title(project_id),
    }

    for proj_id in data:
        results = get_all_files_one_project(proj_id)
        get_files_from_folders(proj_id, results)
        data[proj_id]["files"] = results
    return data


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
parser.add_argument(
    "--size_in_bytes",
    action="store_true",
    help="In TSV output, report file size in bytes instead of MB",
)
options = parser.parse_args()

log = logging.getLogger()
log.setLevel(logging.INFO)

results = get_all_data(options.project_id)

if options.format == "json":
    print(json.dumps(results, indent=2))
elif options.format == "tsv":
    print_tsv(results, size_in_bytes=options.size_in_bytes)
else:
    raise NotImplementedError(f"Uknown output format {options.format}")
