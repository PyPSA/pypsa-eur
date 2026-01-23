# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import copy
from functools import partial, lru_cache

import os, sys, glob
import requests


import pandas as pd
import json

path = workflow.source_path("../scripts/_helpers.py")
sys.path.insert(0, os.path.dirname(path))

from snakemake.utils import update_config


def navigate_config(config, keys, default=None):
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


@lru_cache(maxsize=128)
def _get_config_cached(wildcards_tuple):
    """
    Get full scenario-aware config for given wildcards (internal cached version).

    Parameters
    ----------
    wildcards_tuple : tuple
        Wildcards as frozen tuple for caching

    Returns
    -------
    dict
        Fully resolved config with scenario and wildcard overrides applied
    """
    # Convert back to wildcards dict
    wildcards = dict(wildcards_tuple)

    # Start with base or scenario config
    if config["run"].get("scenarios", {}).get("enable", False) and "run" in wildcards:
        scenario_name = wildcards["run"]
        if scenario_name not in scenarios:
            raise ValueError(
                f"Scenario {scenario_name} not found in {config['run']['scenarios']['file']}"
            )
        base = scenario_config(scenario_name)  # Already cached
    else:
        base = copy.deepcopy(config)

    return base


def get_config(w):
    """
    Get full scenario-aware config for given wildcards.

    This function returns the complete config dictionary with scenario overrides
    (if enabled) and wildcard-based overrides applied.

    Parameters
    ----------
    w : wildcards or dict
        Snakemake wildcards object or dict containing wildcard values

    Returns
    -------
    dict
        Fully resolved config dictionary

    Examples
    --------
    In a Snakemake rule:
        params:
            cfg=lambda w: get_config(w)

    In a helper function:
        def my_function(w):
            cfg = get_config(w)
            return cfg["electricity"]["renewable_carriers"]
    """
    # Convert wildcards to hashable tuple for caching
    wildcards_tuple = tuple(sorted(w.items()))
    return _get_config_cached(wildcards_tuple)


def static_getter(wildcards, keys, default):
    """Getter function for static config values."""
    return navigate_config(config, keys, default)


def dynamic_getter(wildcards, keys, default):
    """Getter function for dynamic config values based on scenario."""
    if "run" not in wildcards.keys():
        return navigate_config(config, keys, default)
    scenario_name = wildcards.run
    if scenario_name not in scenarios:
        raise ValueError(
            f"Scenario {scenario_name} not found in file {config['run']['scenarios']['file']}."
        )
    config_with_scenario = scenario_config(scenario_name)
    return navigate_config(config_with_scenario, keys, default)


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


@lru_cache
def load_data_versions(file_path):
    data_versions = pd.read_csv(
        file_path,
        dtype=str,
        na_filter=False,
        delimiter=",",
        comment="#",
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

    dataset_config = config["data"][
        name
    ]  # TODO as is right now, it is not compatible with config_provider

    # To use PyPSA-Eur as a snakemake module, the path to the versions.csv file needs to be
    # registered relative to the current file with Snakemake:
    fp = workflow.source_path("../data/versions.csv")
    data_versions = load_data_versions(fp)

    dataset = data_versions.loc[
        (data_versions["dataset"] == name)
        & (data_versions["source"] == dataset_config["source"])
        & (data_versions["supported"])  # Limit to supported versions only
        & (
            data_versions["version"] == dataset_config["version"]
            if "latest" != dataset_config["version"]
            else True
        )
        & (data_versions["latest"] if "latest" == dataset_config["version"] else True)
    ]

    if dataset.empty:
        raise ValueError(
            f"Dataset '{name}' with source '{dataset_config['source']}' for '{dataset_config['version']}' not found in data/versions.csv."
        )

    # Return single-row DataFrame as a Series
    dataset = dataset.squeeze()

    # Generate output folder path in the `data` directory
    dataset["folder"] = Path(
        "data", name, dataset["source"], dataset["version"]
    ).as_posix()

    return dataset


def solver_threads(w):
    solver_options = config_provider("solving", "solver_options")(w)
    option_set = config_provider("solving", "solver", "options")(w)
    solver_option_set = solver_options[option_set]
    threads = solver_option_set.get("threads") or solver_option_set.get("Threads") or 4
    return threads


def input_custom_extra_functionality(w):
    path = config_provider(
        "solving", "options", "custom_extra_functionality", default=False
    )(w)
    if path:
        return os.path.join(os.path.dirname(workflow.snakefile), path)
    return []


def solved_previous_horizon(w):
    horizons = config_provider("planning_horizons")(w)
    i = horizons.index(int(w.horizon))
    horizon_p = str(horizons[i - 1])

    return RESULTS + "networks/solved_" + horizon_p + ".nc"


def input_cutout(wildcards, cutout_names="default"):

    cutouts_path = dataset_version("cutout")["folder"]

    if cutout_names == "default":
        cutout_names = config_provider("atlite", "default_cutout")(wildcards)

    if isinstance(cutout_names, list):
        return [f"{cutouts_path}/{cn}.nc" for cn in cutout_names]
    else:
        return f"{cutouts_path}/{cutout_names}.nc"
