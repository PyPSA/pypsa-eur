# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import copy


def get_config(keys, config, default=None):
    """Retrieve a nested value from a dictionary using a tuple of keys."""
    value = config
    for key in keys:
        value = value.get(key, default)
        if value == default:
            return default
    return value


def merge_configs(base_config, scenario_config):
    """Merge base config with a specific scenario without modifying the original."""
    merged = copy.deepcopy(base_config)
    for key, value in scenario_config.items():
        if key in merged and isinstance(merged[key], dict):
            merged[key] = merge_configs(merged[key], value)
        else:
            merged[key] = value
    return merged


def config_provider(*keys, default=None):
    """Dynamically provide config values based on 'run' -> 'name'.

    Usage in Snakemake rules would look something like:
    params:
        my_param=config_provider("key1", "key2", default="some_default_value")
    """

    def static_getter(wildcards):
        """Getter function for static config values."""
        return get_config(keys, config, default)

    def dynamic_getter(wildcards):
        """Getter function for dynamic config values based on scenario."""
        scenario_name = wildcards.run
        if scenario_name not in scenarios:
            raise ValueError(
                f"Scenario {scenario_name} not found in file {config['scenariofile']}."
            )
        merged_config = merge_configs(config, scenarios[scenario_name])
        return get_config(keys, merged_config, default)

    if config["run"].get("scenarios", False):
        return dynamic_getter
    else:
        return static_getter


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
    return f"data/eurostat-energy_balances-june_{report_year}_edition"


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
