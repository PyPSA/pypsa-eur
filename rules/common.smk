# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import copy
from functools import partial, lru_cache

import os, sys, glob
import requests
from tenacity import (
    retry as tenacity_retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)

path = workflow.source_path("../scripts/_helpers.py")
sys.path.insert(0, os.path.dirname(path))

from scripts._helpers import validate_checksum, update_config_from_wildcards
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
def get_full_config(wildcards_tuple):
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

    # Apply wildcard overrides
    result = update_config_from_wildcards(base, wildcards, inplace=False)

    return result


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
    return get_full_config(wildcards_tuple)


def static_getter(wildcards, keys, default):
    """Getter function for static config values."""
    config_with_wildcards = update_config_from_wildcards(
        config, wildcards, inplace=False
    )
    return navigate_config(config_with_wildcards, keys, default)


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
    config_with_wildcards = update_config_from_wildcards(
        config_with_scenario, wildcards, inplace=False
    )
    return navigate_config(config_with_wildcards, keys, default)


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


@tenacity_retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=4, max=10),
    retry=retry_if_exception_type(
        (requests.HTTPError, requests.ConnectionError, requests.Timeout)
    ),
)
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
    # Send a HEAD request to avoid fetching full response
    response = requests.head(url, timeout=timeout, allow_redirects=True)
    # Raise HTTPError for transient errors
    # 429: Too Many Requests (rate limiting)
    # 500, 502, 503, 504: Server errors
    if response.status_code in (429, 500, 502, 503, 504):
        response.raise_for_status()
    return response.status_code == 200


def solved_previous_horizon(w):
    horizons = config_provider("planning_horizons")(w)
    i = horizons.index(int(w.horizon))
    planning_horizon_p = str(horizons[i - 1])

    return RESULTS + "networks/solved_" + planning_horizon_p + ".nc"


def input_cutout(wildcards, cutout_names="default"):
    if cutout_names == "default":
        cutout_names = config_provider("atlite", "default_cutout")(wildcards)
    if isinstance(cutout_names, list):
        return [CDIR.joinpath(cn + ".nc").as_posix() for cn in cutout_names]
    else:
        return CDIR.joinpath(cutout_names + ".nc").as_posix()
