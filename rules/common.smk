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
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


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
