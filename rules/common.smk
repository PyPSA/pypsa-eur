# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import os, sys, glob

helper_source_path = [match for match in glob.glob("**/_helpers.py", recursive=True)]

for path in helper_source_path:
    path = os.path.dirname(os.path.abspath(path))
    sys.path.insert(0, os.path.abspath(path))

from _helpers import validate_checksum


def solver_threads(w):
    solver_options = config["solving"]["solver_options"]
    option_set = config["solving"]["solver"]["options"]
    threads = solver_options[option_set].get("threads", 4)
    return threads


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m") or w.clusters.endswith("c"):
        return int(factor * (55000 + 600 * int(w.clusters[:-1])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


def input_custom_extra_functionality(w):
    path = config["solving"]["options"].get("custom_extra_functionality", False)
    if path:
        return os.path.join(os.path.dirname(workflow.snakefile), path)
    return []


# Check if the workflow has access to the internet by trying to access the HEAD of specified url
def has_internet_access(url="www.zenodo.org") -> bool:
    import http.client as http_client

    # based on answer and comments from
    # https://stackoverflow.com/a/29854274/11318472
    conn = http_client.HTTPConnection(url, timeout=5)  # need access to zenodo anyway
    try:
        conn.request("HEAD", "/")
        return True
    except:
        return False
    finally:
        conn.close()


def input_eurostat(w):
    # 2016 includes BA, 2017 does not
    report_year = config["energy"]["eurostat_report_year"]
    return f"data/bundle-sector/eurostat-energy_balances-june_{report_year}_edition"


def solved_previous_horizon(wildcards):
    planning_horizons = config["scenario"]["planning_horizons"]
    i = planning_horizons.index(int(wildcards.planning_horizons))
    planning_horizon_p = str(planning_horizons[i - 1])
    return (
        RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_"
        + planning_horizon_p
        + ".nc"
    )
