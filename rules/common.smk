# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


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
        + "postnetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_"
        + planning_horizon_p
        + ".nc"
    )
