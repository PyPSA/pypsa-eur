# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import contextlib
import copy
import hashlib
import logging
import os
import re
import urllib
from functools import partial
from os.path import exists
from pathlib import Path
from shutil import copyfile

import pandas as pd
import pytz
import requests
import yaml
from snakemake.utils import update_config
from tqdm import tqdm

logger = logging.getLogger(__name__)

REGION_COLS = ["geometry", "name", "x", "y", "country"]


def copy_default_files(workflow):
    default_files = {
        "config/config.default.yaml": "config/config.yaml",
        "config/scenarios.template.yaml": "config/scenarios.yaml",
    }
    for template, target in default_files.items():
        target = os.path.join(workflow.current_basedir, target)
        template = os.path.join(workflow.current_basedir, template)
        if not exists(target) and exists(template):
            copyfile(template, target)


def get_scenarios(run):
    scenario_config = run.get("scenarios", {})
    if run["name"] and scenario_config.get("enable"):
        fn = Path(scenario_config["file"])
        if fn.exists():
            scenarios = yaml.safe_load(fn.read_text())
            if run["name"] == "all":
                run["name"] = list(scenarios.keys())
            return scenarios
    return {}


def get_rdir(run):
    scenario_config = run.get("scenarios", {})
    if run["name"] and scenario_config.get("enable"):
        RDIR = "{run}/"
    elif run["name"]:
        RDIR = run["name"] + "/"
    else:
        RDIR = ""

    prefix = run.get("prefix", "")
    if prefix:
        RDIR = f"{prefix}/{RDIR}"

    return RDIR


def get_run_path(fn, dir, rdir, shared_resources, exclude_from_shared):
    """
    Dynamically provide paths based on shared resources and filename.

    Use this function for snakemake rule inputs or outputs that should be
    optionally shared across runs or created individually for each run.

    Parameters
    ----------
    fn : str
        The filename for the path to be generated.
    dir : str
        The base directory.
    rdir : str
        Relative directory for non-shared resources.
    shared_resources : str or bool
        Specifies which resources should be shared.
        - If string is "base", special handling for shared "base" resources (see notes).
        - If random string other than "base", this folder is used instead of the `rdir` keyword.
        - If boolean, directly specifies if the resource is shared.
    exclude_from_shared: list
        List of filenames to exclude from shared resources. Only relevant if shared_resources is "base".

    Returns
    -------
    str
        Full path where the resource should be stored.

    Notes
    -----
    Special case for "base" allows no wildcards other than "technology", "year"
    and "scope" and excludes filenames starting with "networks/elec" or
    "add_electricity". All other resources are shared.
    """
    if shared_resources == "base":
        pattern = r"\{([^{}]+)\}"
        existing_wildcards = set(re.findall(pattern, fn))
        irrelevant_wildcards = {"technology", "year", "scope", "kind"}
        no_relevant_wildcards = not existing_wildcards - irrelevant_wildcards
        not_shared_rule = (
            not fn.startswith("networks/elec")
            and not fn.startswith("add_electricity")
            and not any(fn.startswith(ex) for ex in exclude_from_shared)
        )
        is_shared = no_relevant_wildcards and not_shared_rule
        rdir = "" if is_shared else rdir
    elif isinstance(shared_resources, str):
        rdir = shared_resources + "/"
    elif isinstance(shared_resources, bool):
        rdir = "" if shared_resources else rdir
    else:
        raise ValueError(
            "shared_resources must be a boolean, str, or 'base' for special handling."
        )

    return f"{dir}{rdir}{fn}"


def path_provider(dir, rdir, shared_resources, exclude_from_shared):
    """
    Returns a partial function that dynamically provides paths based on shared
    resources and the filename.

    Returns
    -------
    partial function
        A partial function that takes a filename as input and
        returns the path to the file based on the shared_resources parameter.
    """
    return partial(
        get_run_path,
        dir=dir,
        rdir=rdir,
        shared_resources=shared_resources,
        exclude_from_shared=exclude_from_shared,
    )


def get_opt(opts, expr, flags=None):
    """
    Return the first option matching the regular expression.

    The regular expression is case-insensitive by default.
    """
    if flags is None:
        flags = re.IGNORECASE
    for o in opts:
        match = re.match(expr, o, flags=flags)
        if match:
            return match.group(0)
    return None


def find_opt(opts, expr):
    """
    Return if available the float after the expression.
    """
    for o in opts:
        if expr in o:
            m = re.findall(r"m?\d+(?:[\.p]\d+)?", o)
            if len(m) > 0:
                return True, float(m[-1].replace("p", ".").replace("m", "-"))
            else:
                return True, None
    return False, None


# Define a context manager to temporarily mute print statements
@contextlib.contextmanager
def mute_print():
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull):
            yield


def set_scenario_config(snakemake):
    scenario = snakemake.config["run"].get("scenarios", {})
    if scenario.get("enable") and "run" in snakemake.wildcards.keys():
        try:
            with open(scenario["file"], "r") as f:
                scenario_config = yaml.safe_load(f)
        except FileNotFoundError:
            # fallback for mock_snakemake
            script_dir = Path(__file__).parent.resolve()
            root_dir = script_dir.parent
            with open(root_dir / scenario["file"], "r") as f:
                scenario_config = yaml.safe_load(f)
        update_config(snakemake.config, scenario_config[snakemake.wildcards.run])


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """
    import logging
    import sys

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)

    # Setup a function to handle uncaught exceptions and include them with their stacktrace into logfiles
    def handle_exception(exc_type, exc_value, exc_traceback):
        # Log the exception
        logger = logging.getLogger()
        logger.error(
            "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
        )

    sys.excepthook = handle_exception


def update_p_nom_max(n):
    # if extendable carriers (solar/onwind/...) have capacity >= 0,
    # e.g. existing assets from the OPSD project are included to the network,
    # the installed capacity might exceed the expansion limit.
    # Hence, we update the assumptions.

    n.generators.p_nom_max = n.generators[["p_nom_min", "p_nom_max"]].max(1)


def aggregate_p_nom(n):
    return pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum(),
            n.storage_units.groupby("carrier").p_nom_opt.sum(),
            n.links.groupby("carrier").p_nom_opt.sum(),
            n.loads_t.p.groupby(n.loads.carrier, axis=1).sum().mean(),
        ]
    )


def aggregate_p(n):
    return pd.concat(
        [
            n.generators_t.p.sum().groupby(n.generators.carrier).sum(),
            n.storage_units_t.p.sum().groupby(n.storage_units.carrier).sum(),
            n.stores_t.p.sum().groupby(n.stores.carrier).sum(),
            -n.loads_t.p.sum().groupby(n.loads.carrier).sum(),
        ]
    )


def get(item, investment_year=None):
    """
    Check whether item depends on investment year.
    """
    if not isinstance(item, dict):
        return item
    elif investment_year in item.keys():
        return item[investment_year]
    else:
        logger.warning(
            f"Investment key {investment_year} not found in dictionary {item}."
        )
        keys = sorted(item.keys())
        if investment_year < keys[0]:
            logger.warning(f"Lower than minimum key. Taking minimum key {keys[0]}")
            return item[keys[0]]
        elif investment_year > keys[-1]:
            logger.warning(f"Higher than maximum key. Taking maximum key {keys[0]}")
            return item[keys[-1]]
        else:
            logger.warning(
                "Interpolate linearly between the next lower and next higher year."
            )
            lower_key = max(k for k in keys if k < investment_year)
            higher_key = min(k for k in keys if k > investment_year)
            lower = item[lower_key]
            higher = item[higher_key]
            return lower + (higher - lower) * (investment_year - lower_key) / (
                higher_key - lower_key
            )


def aggregate_e_nom(n):
    return pd.concat(
        [
            (n.storage_units["p_nom_opt"] * n.storage_units["max_hours"])
            .groupby(n.storage_units["carrier"])
            .sum(),
            n.stores["e_nom_opt"].groupby(n.stores.carrier).sum(),
        ]
    )


def aggregate_p_curtailed(n):
    return pd.concat(
        [
            (
                (
                    n.generators_t.p_max_pu.sum().multiply(n.generators.p_nom_opt)
                    - n.generators_t.p.sum()
                )
                .groupby(n.generators.carrier)
                .sum()
            ),
            (
                (n.storage_units_t.inflow.sum() - n.storage_units_t.p.sum())
                .groupby(n.storage_units.carrier)
                .sum()
            ),
        ]
    )


def aggregate_costs(n, flatten=False, opts=None, existing_only=False):
    components = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
        Transformer=("s_nom", None),
    )

    costs = {}
    for c, (p_nom, p_attr) in zip(
        n.iterate_components(components.keys(), skip_empty=False), components.values()
    ):
        if c.df.empty:
            continue
        if not existing_only:
            p_nom += "_opt"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs)

    if flatten:
        assert opts is not None
        conv_techs = opts["conv_techs"]

        costs = costs.reset_index(level=0, drop=True)
        costs = costs["capital"].add(
            costs["marginal"].rename({t: t + " marginal" for t in conv_techs}),
            fill_value=0.0,
        )

    return costs


def progress_retrieve(url, file, disable=False):
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    # Hotfix - Bug, tqdm not working with disable=False
    disable = True

    if disable:
        response = requests.get(url, headers=headers, stream=True)
        with open(file, "wb") as f:
            f.write(response.content)
    else:
        response = requests.get(url, headers=headers, stream=True)
        total_size = int(response.headers.get("content-length", 0))
        chunk_size = 1024

        with tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=str(file),
        ) as t:
            with open(file, "wb") as f:
                for data in response.iter_content(chunk_size=chunk_size):
                    f.write(data)
                    t.update(len(data))


def mock_snakemake(
    rulename,
    root_dir=None,
    configfiles=None,
    submodule_dir="workflow/submodules/pypsa-eur",
    **wildcards,
):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    root_dir: str/path-like
        path to the root directory of the snakemake project
    configfiles: list, str
        list of configfiles to be used to update the config
    submodule_dir: str, Path
        in case PyPSA-Eur is used as a submodule, submodule_dir is
        the path of pypsa-eur relative to the project directory.
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm
    from pypsa.definitions.structures import Dict
    from snakemake.api import Workflow
    from snakemake.common import SNAKEFILE_CHOICES
    from snakemake.script import Snakemake
    from snakemake.settings.types import (
        ConfigSettings,
        DAGSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

    script_dir = Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = Path(root_dir).resolve()

    user_in_script_dir = Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif Path.cwd().resolve() != root_dir:
        raise RuntimeError(
            "mock_snakemake has to be run from the repository root"
            f" {root_dir} or scripts directory {script_dir}"
        )
    try:
        for p in SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break
        if configfiles is None:
            configfiles = []
        elif isinstance(configfiles, str):
            configfiles = [configfiles]

        resource_settings = ResourceSettings()
        config_settings = ConfigSettings(configfiles=map(Path, configfiles))
        workflow_settings = WorkflowSettings()
        storage_settings = StorageSettings()
        dag_settings = DAGSettings(rerun_triggers=[])
        workflow = Workflow(
            config_settings,
            resource_settings,
            workflow_settings,
            storage_settings,
            dag_settings,
            storage_provider_settings=dict(),
        )
        workflow.include(snakefile)

        if configfiles:
            for f in configfiles:
                if not os.path.exists(f):
                    raise FileNotFoundError(f"Config file {f} does not exist.")
                workflow.configfile(f)

        workflow.global_resources = {}
        rule = workflow.get_rule(rulename)
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i, _ in enumerate(io):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    finally:
        if user_in_script_dir:
            os.chdir(script_dir)
    return snakemake


def generate_periodic_profiles(dt_index, nodes, weekly_profile, localize=None):
    """
    Give a 24*7 long list of weekly hourly profiles, generate this for each
    country for the period dt_index, taking account of time zones and summer
    time.
    """
    weekly_profile = pd.Series(weekly_profile, range(24 * 7))

    week_df = pd.DataFrame(index=dt_index, columns=nodes)

    for node in nodes:
        ct = node[:2] if node[:2] != "XK" else "RS"
        timezone = pytz.timezone(pytz.country_timezones[ct][0])
        tz_dt_index = dt_index.tz_convert(timezone)
        week_df[node] = [24 * dt.weekday() + dt.hour for dt in tz_dt_index]
        week_df[node] = week_df[node].map(weekly_profile)

    week_df = week_df.tz_localize(localize)

    return week_df


def parse(infix):
    """
    Recursively parse a chained wildcard expression into a dictionary or a YAML
    object.

    Parameters
    ----------
    list_to_parse : list
        The list to parse.

    Returns
    -------
    dict or YAML object
        The parsed list.
    """
    if len(infix) == 1:
        return yaml.safe_load(infix[0])
    else:
        return {infix.pop(0): parse(infix)}


def update_config_from_wildcards(config, w, inplace=True):
    """
    Parses configuration settings from wildcards and updates the config.
    """

    if not inplace:
        config = copy.deepcopy(config)

    if w.get("opts"):
        opts = w.opts.split("-")

        if nhours := get_opt(opts, r"^\d+(h|seg)$"):
            config["clustering"]["temporal"]["resolution_elec"] = nhours

        co2l_enable, co2l_value = find_opt(opts, "Co2L")
        if co2l_enable:
            config["electricity"]["co2limit_enable"] = True
            if co2l_value is not None:
                config["electricity"]["co2limit"] = (
                    co2l_value * config["electricity"]["co2base"]
                )

        gasl_enable, gasl_value = find_opt(opts, "CH4L")
        if gasl_enable:
            config["electricity"]["gaslimit_enable"] = True
            if gasl_value is not None:
                config["electricity"]["gaslimit"] = gasl_value * 1e6

        if "Ept" in opts:
            config["costs"]["emission_prices"]["co2_monthly_prices"] = True

        ep_enable, ep_value = find_opt(opts, "Ep")
        if ep_enable:
            config["costs"]["emission_prices"]["enable"] = True
            if ep_value is not None:
                config["costs"]["emission_prices"]["co2"] = ep_value

        if "ATK" in opts:
            config["autarky"]["enable"] = True
            if "ATKc" in opts:
                config["autarky"]["by_country"] = True

        attr_lookup = {
            "p": "p_nom_max",
            "e": "e_nom_max",
            "c": "capital_cost",
            "m": "marginal_cost",
        }
        for o in opts:
            flags = ["+e", "+p", "+m", "+c"]
            if all(flag not in o for flag in flags):
                continue
            carrier, attr_factor = o.split("+")
            attr = attr_lookup[attr_factor[0]]
            factor = float(attr_factor[1:])
            if not isinstance(config["adjustments"]["electricity"], dict):
                config["adjustments"]["electricity"] = dict()
            update_config(
                config["adjustments"]["electricity"], {attr: {carrier: factor}}
            )

    if w.get("sector_opts"):
        opts = w.sector_opts.split("-")

        if "T" in opts:
            config["sector"]["transport"] = True

        if "H" in opts:
            config["sector"]["heating"] = True

        if "B" in opts:
            config["sector"]["biomass"] = True

        if "I" in opts:
            config["sector"]["industry"] = True

        if "A" in opts:
            config["sector"]["agriculture"] = True

        if "CCL" in opts:
            config["solving"]["constraints"]["CCL"] = True

        eq_value = get_opt(opts, r"^EQ+\d*\.?\d+(c|)")
        for o in opts:
            if eq_value is not None:
                config["solving"]["constraints"]["EQ"] = eq_value
            elif "EQ" in o:
                config["solving"]["constraints"]["EQ"] = True
            break

        if "BAU" in opts:
            config["solving"]["constraints"]["BAU"] = True

        if "SAFE" in opts:
            config["solving"]["constraints"]["SAFE"] = True

        if nhours := get_opt(opts, r"^\d+(h|sn|seg)$"):
            config["clustering"]["temporal"]["resolution_sector"] = nhours

        if "decentral" in opts:
            config["sector"]["electricity_transmission_grid"] = False

        if "noH2network" in opts:
            config["sector"]["H2_network"] = False

        if "nowasteheat" in opts:
            config["sector"]["use_fischer_tropsch_waste_heat"] = False
            config["sector"]["use_methanolisation_waste_heat"] = False
            config["sector"]["use_haber_bosch_waste_heat"] = False
            config["sector"]["use_methanation_waste_heat"] = False
            config["sector"]["use_fuel_cell_waste_heat"] = False
            config["sector"]["use_electrolysis_waste_heat"] = False

        if "nodistrict" in opts:
            config["sector"]["district_heating"]["progress"] = 0.0

        dg_enable, dg_factor = find_opt(opts, "dist")
        if dg_enable:
            config["sector"]["electricity_distribution_grid"] = True
            if dg_factor is not None:
                config["sector"][
                    "electricity_distribution_grid_cost_factor"
                ] = dg_factor

        if "biomasstransport" in opts:
            config["sector"]["biomass_transport"] = True

        _, maxext = find_opt(opts, "linemaxext")
        if maxext is not None:
            config["lines"]["max_extension"] = maxext * 1e3
            config["links"]["max_extension"] = maxext * 1e3

        _, co2l_value = find_opt(opts, "Co2L")
        if co2l_value is not None:
            config["co2_budget"] = float(co2l_value)

        if co2_distribution := get_opt(opts, r"^(cb)\d+(\.\d+)?(ex|be)$"):
            config["co2_budget"] = co2_distribution

        if co2_budget := get_opt(opts, r"^(cb)\d+(\.\d+)?$"):
            config["co2_budget"] = float(co2_budget[2:])

        attr_lookup = {
            "p": "p_nom_max",
            "e": "e_nom_max",
            "c": "capital_cost",
            "m": "marginal_cost",
        }
        for o in opts:
            flags = ["+e", "+p", "+m", "+c"]
            if all(flag not in o for flag in flags):
                continue
            carrier, attr_factor = o.split("+")
            attr = attr_lookup[attr_factor[0]]
            factor = float(attr_factor[1:])
            if not isinstance(config["adjustments"]["sector"], dict):
                config["adjustments"]["sector"] = dict()
            update_config(config["adjustments"]["sector"], {attr: {carrier: factor}})

        _, sdr_value = find_opt(opts, "sdr")
        if sdr_value is not None:
            config["costs"]["social_discountrate"] = sdr_value / 100

        _, seq_limit = find_opt(opts, "seq")
        if seq_limit is not None:
            config["sector"]["co2_sequestration_potential"] = seq_limit

        # any config option can be represented in wildcard
        for o in opts:
            if o.startswith("CF+"):
                infix = o.split("+")[1:]
                update_config(config, parse(infix))

    if not inplace:
        return config


def get_checksum_from_zenodo(file_url):
    parts = file_url.split("/")
    record_id = parts[parts.index("records") + 1]
    filename = parts[-1]

    response = requests.get(f"https://zenodo.org/api/records/{record_id}", timeout=30)
    response.raise_for_status()
    data = response.json()

    for file in data["files"]:
        if file["key"] == filename:
            return file["checksum"]
    return None


def validate_checksum(file_path, zenodo_url=None, checksum=None):
    """
    Validate file checksum against provided or Zenodo-retrieved checksum.
    Calculates the hash of a file using 64KB chunks. Compares it against a
    given checksum or one from a Zenodo URL.

    Parameters
    ----------
    file_path : str
        Path to the file for checksum validation.
    zenodo_url : str, optional
        URL of the file on Zenodo to fetch the checksum.
    checksum : str, optional
        Checksum (format 'hash_type:checksum_value') for validation.

    Raises
    ------
    AssertionError
        If the checksum does not match, or if neither `checksum` nor `zenodo_url` is provided.


    Examples
    --------
    >>> validate_checksum("/path/to/file", checksum="md5:abc123...")
    >>> validate_checksum(
    ...     "/path/to/file",
    ...     zenodo_url="https://zenodo.org/records/12345/files/example.txt",
    ... )

    If the checksum is invalid, an AssertionError will be raised.
    """
    assert checksum or zenodo_url, "Either checksum or zenodo_url must be provided"
    if zenodo_url:
        checksum = get_checksum_from_zenodo(zenodo_url)
    hash_type, checksum = checksum.split(":")
    hasher = hashlib.new(hash_type)
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):  # 64kb chunks
            hasher.update(chunk)
    calculated_checksum = hasher.hexdigest()
    assert (
        calculated_checksum == checksum
    ), "Checksum is invalid. This may be due to an incomplete download. Delete the file and re-execute the rule."


def get_snapshots(snapshots, drop_leap_day=False, freq="h", **kwargs):
    """
    Returns pandas DateTimeIndex potentially without leap days.
    """

    time = pd.date_range(freq=freq, **snapshots, **kwargs)
    if drop_leap_day and time.is_leap_year.any():
        time = time[~((time.month == 2) & (time.day == 29))]

    return time


def rename_techs(label: str) -> str:
    """
    Rename technology labels for better readability.

    Removes some prefixes and renames if certain conditions defined in function body are met.

    Parameters:
    ----------
    label: str
        Technology label to be renamed

    Returns:
    -------
    str
        Renamed label
    """
    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral ",
    ]

    rename_if_contains = [
        "CHP",
        "gas boiler",
        "biogas",
        "solar thermal",
        "air heat pump",
        "ground heat pump",
        "resistive heater",
        "Fischer-Tropsch",
    ]

    rename_if_contains_dict = {
        "water tanks": "hot water storage",
        "retrofitting": "building retrofitting",
        # "H2 Electrolysis": "hydrogen storage",
        # "H2 Fuel Cell": "hydrogen storage",
        # "H2 pipeline": "hydrogen storage",
        "battery": "battery storage",
        "H2 for industry": "H2 for industry",
        "land transport fuel cell": "land transport fuel cell",
        "land transport oil": "land transport oil",
        "oil shipping": "shipping oil",
        # "CC": "CC"
    }

    rename = {
        "solar": "solar PV",
        "Sabatier": "methanation",
        "offwind": "offshore wind",
        "offwind-ac": "offshore wind (AC)",
        "offwind-dc": "offshore wind (DC)",
        "offwind-float": "offshore wind (Float)",
        "onwind": "onshore wind",
        "ror": "hydroelectricity",
        "hydro": "hydroelectricity",
        "PHS": "hydroelectricity",
        "NH3": "ammonia",
        "co2 Store": "DAC",
        "co2 stored": "CO2 sequestration",
        "AC": "transmission lines",
        "DC": "transmission lines",
        "B2B": "transmission lines",
    }

    for ptr in prefix_to_remove:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old, new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old, new in rename.items():
        if old == label:
            label = new
    return label
