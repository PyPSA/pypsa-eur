# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import copy
from functools import partial, lru_cache

import os, sys, glob

import pandas as pd
import json

path = workflow.source_path("../scripts/_helpers.py")
sys.path.insert(0, os.path.dirname(path))

from scripts._helpers import validate_checksum, update_config_from_wildcards
from snakemake.utils import update_config


def get_config(config, keys, default=None):
    """Retrieve a nested value from a dictionary using a tuple of keys."""
    value = config
    for key in keys:
        if isinstance(value, list):
            value = value[key]
        else:
            value = value.get(key, default)
        if value == default:
            return default
    return value


def merge_configs(base_config, scenario_config):
    """Merge base config with a specific scenario without modifying the original."""
    merged = copy.deepcopy(base_config)
    update_config(merged, scenario_config)
    return merged


@lru_cache
def scenario_config(scenario_name):
    """Retrieve a scenario config based on the overrides from the scenario file."""
    return merge_configs(config, scenarios[scenario_name])


def static_getter(wildcards, keys, default):
    """Getter function for static config values."""
    config_with_wildcards = update_config_from_wildcards(
        config, wildcards, inplace=False
    )
    return get_config(config_with_wildcards, keys, default)


def dynamic_getter(wildcards, keys, default):
    """Getter function for dynamic config values based on scenario."""
    if "run" not in wildcards.keys():
        return get_config(config, keys, default)
    scenario_name = wildcards.run
    if scenario_name not in scenarios:
        raise ValueError(
            f"Scenario {scenario_name} not found in file {config['run']['scenarios']['file']}."
        )
    config_with_scenario = scenario_config(scenario_name)
    config_with_wildcards = update_config_from_wildcards(
        config_with_scenario, wildcards, inplace=False
    )
    return get_config(config_with_wildcards, keys, default)


def config_provider(*keys, default=None):
    """Dynamically provide config values based on 'run' -> 'name'.

    Usage in Snakemake rules would look something like:
    params:
        my_param=config_provider("key1", "key2", default="some_default_value")
    """
    # Using functools.partial to freeze certain arguments in our getter functions.
    if config["run"].get("scenarios", {}).get("enable", False):
        return partial(dynamic_getter, keys=keys, default=default)
    else:
        return partial(static_getter, keys=keys, default=default)


def dataset_version(
    name: str,
) -> pd.Series:
    """
    Return the dataset version information and url for a given dataset name.

    The dataset name is used to determine the source and version of the dataset from the configuration.
    Then the 'data/versions.csv' file is queried to find the matching dataset entry.

    Parameters:
    name: str
        The name of the dataset to retrieve version information for.

    Returns:
    pd.Series
        A pandas Series containing the dataset version information, including source, version, tags, and URL
    """

    @lru_cache
    def load_data_versions(file_path):
        data_versions = pd.read_csv(
            file_path, dtype=str, na_filter=False, delimiter=",", comment="#"
        )

        # Turn 'tags' column from string representation of list to individual columns
        data_versions["tags"] = data_versions["tags"].apply(
            lambda x: json.loads(x.replace("'", '"'))
        )
        exploded = data_versions.explode("tags")
        dummies = pd.get_dummies(exploded["tags"], dtype=bool)
        tags_matrix = dummies.groupby(dummies.index).max()
        data_versions = data_versions.join(tags_matrix)

        return data_versions

    dataset_config = config["data"][
        name
    ]  # TODO as is right now, it is not compatible with config_provider
    data_versions = load_data_versions("data/versions.csv")

    dataset = data_versions.loc[
        (data_versions["dataset"] == name)
        & (data_versions["source"] == dataset_config["source"])
        & (data_versions["supported"])  # Limit to supported versions only
        & (
            data_versions["version"] == dataset_config["version_or_latest"]
            if "latest" != dataset_config["version_or_latest"]
            else True
        )
        & (
            data_versions["latest"]
            if "latest" == dataset_config["version_or_latest"]
            else True
        )
    ]

    if dataset.empty:
        raise ValueError(
            f"Dataset '{name}' with source '{dataset_config['source']}' for '{dataset_config['version_or_latest']}' not found in data/versions.csv."
        )

    # Return single-row DataFrame as a Series
    dataset = dataset.squeeze()

    # Generate output folder path in the `data` directory
    dataset["folder"] = Path("data", name, dataset["source"], dataset["version"])

    return dataset


def handle_data_requests(params, output):
    import os
    import requests
    from zipfile import ZipFile
    from pathlib import Path

    response = requests.get(params["url"])
    with open(output.zip, "wb") as f:
        f.write(response.content)

    output_folder = Path(output["zip"]).parent
    unpack_archive(output.zip, output_folder)


def solver_threads(w):
    solver_options = config_provider("solving", "solver_options")(w)
    option_set = config_provider("solving", "solver", "options")(w)
    solver_option_set = solver_options[option_set]
    threads = solver_option_set.get("threads") or solver_option_set.get("Threads") or 4
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
    if w.clusters == "all" or w.clusters == "adm":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


def input_custom_extra_functionality(w):
    path = config_provider(
        "solving", "options", "custom_extra_functionality", default=False
    )(w)
    if path:
        return os.path.join(os.path.dirname(workflow.snakefile), path)
    return []


def has_internet_access(url: str = "https://www.zenodo.org", timeout: int = 5) -> bool:
    """
    Checks if internet connection is available by sending a HEAD request
    to a reliable server like Zenodo.

    Parameters:
    - url (str): The URL to check for internet connection. Default is Zenodo.
    - timeout (int | float): The maximum time (in seconds) the request should wait.

    Returns:
    - bool: True if the internet is available, otherwise False.
    """
    try:
        # Send a HEAD request to avoid fetching full response
        response = requests.head(url, timeout=timeout, allow_redirects=True)
        return response.status_code == 200
    except requests.ConnectionError:  # (e.g., no internet, DNS issues)
        return False
    except requests.Timeout:  # (e.g., slow or no network)
        return False


def solved_previous_horizon(w):
    planning_horizons = config_provider("scenario", "planning_horizons")(w)
    i = planning_horizons.index(int(w.planning_horizons))
    planning_horizon_p = str(planning_horizons[i - 1])

    return (
        RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_"
        + planning_horizon_p
        + ".nc"
    )


def input_cutout(wildcards, cutout_names="default"):

    CDIR = dataset_version("cutout")["folder"]

    if cutout_names == "default":
        cutout_names = config_provider("atlite", "default_cutout")(wildcards)

    if isinstance(cutout_names, list):
        return [(CDIR / f"{cn}.nc").as_posix() for cn in cutout_names]
    else:
        return (CDIR / f"{cutout_names}.nc").as_posix()
