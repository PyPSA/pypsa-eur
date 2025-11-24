# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


"""
Prepare PyPSA network for solving with various operational constraints and
temporal adjustments.

- adding an annual **limit** of carbon-dioxide emissions,
- adding an exogenous **price** per tonne emissions of carbon-dioxide (or other kinds),
- setting an **N-1 security margin** factor for transmission line capacities,
- specifying an expansion limit on the **cost** of transmission expansion,
- specifying an expansion limit on the **volume** of transmission expansion, and
- reducing the **temporal** resolution by averaging over multiple hours
  or segmenting time series into chunks of varying lengths using ``tsam``.

Description
-----------

.. note::

    This script's functionality has been integrated into :mod:`compose_network`
    in the streamlined workflow. This script is maintained for backwards compatibility.

.. tip::
    The deprecated rule :mod:`prepare_elec_networks` has been replaced by the
    unified :mod:`compose_network` script which handles network preparation for
    all scenarios.
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import (
    PYPSA_V1,
    configure_logging,
    get,
    load_costs,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import set_transmission_costs

# Allow for PyPSA versions <0.35
if PYPSA_V1:
    from pypsa.common import expand_series
else:
    from pypsa.descriptors import expand_series


idx = pd.IndexSlice

logger = logging.getLogger(__name__)

# Conversion constant for CO2 emissions
GT_TO_TONNES = 1e9  # Gigatonnes to tonnes conversion


def modify_attribute(n, adjustments, investment_year, modification="factor"):
    if not adjustments[modification]:
        return
    change_dict = adjustments[modification]
    for c in change_dict.keys():
        if c not in n.component_attrs.keys():
            logger.warning(f"{c} needs to be a PyPSA Component")
            continue
        for carrier in change_dict[c].keys():
            ind_i = n.df(c)[n.df(c).carrier == carrier].index
            if ind_i.empty:
                continue
            for parameter in change_dict[c][carrier].keys():
                if parameter not in n.df(c).columns:
                    logger.warning(f"Attribute {parameter} needs to be in {c} columns.")
                    continue
                if investment_year:
                    factor = get(change_dict[c][carrier][parameter], investment_year)
                else:
                    factor = change_dict[c][carrier][parameter]
                if modification == "factor":
                    logger.info(f"Modify {parameter} of {carrier} by factor {factor} ")
                    n.df(c).loc[ind_i, parameter] *= factor
                elif modification == "absolute":
                    logger.info(f"Set {parameter} of {carrier} to {factor} ")
                    n.df(c).loc[ind_i, parameter] = factor
                else:
                    logger.warning(
                        f"{modification} needs to be either 'absolute' or 'factor'."
                    )


def maybe_adjust_costs_and_potentials(n, adjustments, investment_year=None):
    if not adjustments:
        return
    for modification in adjustments.keys():
        modify_attribute(n, adjustments, investment_year, modification)


def add_co2limit(
    n: pypsa.Network,
    co2_max: float,
    co2_min: float | None = None,
    nyears: float = 1.0,
    suffix: str = "",
) -> None:
    """
    Add global CO2 emissions constraint(s) to the network.

    Parameters
    ----------
    n : pypsa.Network
        Network to add constraints to
    co2_max : float
        Maximum CO2 emissions in Gt CO2. For perfect foresight with investment_period:
        yearly limit for that period. For perfect foresight without investment_period:
        total budget over entire planning horizon. For overnight/myopic: annual limit per year.
    co2_min : float or None, default None
        Minimum CO2 emissions in Gt CO2. For perfect foresight with investment_period:
        yearly limit for that period. For perfect foresight without investment_period:
        total budget over entire planning horizon. For overnight/myopic: annual limit per year.
    nyears : float, default 1.0
        Number of years for overnight/myopic optimization (ignored for perfect foresight)
    suffix : str, default ""
        Suffix to add to constraint names

    """
    # FAIL FAST: Validate inputs
    if pd.isna(co2_max):
        raise ValueError(
            f"co2_max cannot be NaN. Received: {co2_max}. "
            "Ensure CO2 upper constraint is properly configured."
        )

    if co2_min is not None and pd.isna(co2_min):
        raise ValueError(
            f"co2_min cannot be NaN. Received: {co2_min}. "
            "Ensure CO2 lower constraint is properly configured."
        )

    n.add(
        "GlobalConstraint",
        "CO2Limit" + suffix,
        type="co2_limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=co2_max * GT_TO_TONNES * nyears,
    )

    if co2_min is not None:
        n.add(
            "GlobalConstraint",
            "CO2Min" + suffix,
            type="co2_limit",
            carrier_attribute="co2_emissions",
            sense=">=",
            constant=co2_min * GT_TO_TONNES * nyears,
        )


def add_gaslimit(n, gaslimit, Nyears=1.0):
    sel = n.carriers.index.intersection(["OCGT", "CCGT", "CHP"])
    n.carriers.loc[sel, "gas_usage"] = 1.0

    n.add(
        "GlobalConstraint",
        "GasLimit",
        carrier_attribute="gas_usage",
        sense="<=",
        constant=gaslimit * Nyears,
    )


def add_emission_prices(n, emission_prices={"co2": 0.0}, exclude_co2=False):
    if exclude_co2:
        emission_prices.pop("co2")
    ep = (
        pd.Series(emission_prices).rename(lambda x: x + "_emissions")
        * n.carriers.filter(like="_emissions")
    ).sum(axis=1)
    gen_ep = n.generators.carrier.map(ep) / n.generators.efficiency
    n.generators["marginal_cost"] += gen_ep
    n.generators_t["marginal_cost"] += gen_ep[n.generators_t["marginal_cost"].columns]
    su_ep = n.storage_units.carrier.map(ep) / n.storage_units.efficiency_dispatch
    n.storage_units["marginal_cost"] += su_ep


def add_dynamic_emission_prices(n, fn):
    co2_price = pd.read_csv(fn, index_col=0, parse_dates=True)
    co2_price = co2_price[~co2_price.index.duplicated()]
    co2_price = co2_price.reindex(n.snapshots).ffill().bfill()

    emissions = (
        n.generators.carrier.map(n.carriers.co2_emissions) / n.generators.efficiency
    )
    co2_cost = expand_series(emissions, n.snapshots).T.mul(co2_price.iloc[:, 0], axis=0)

    static = n.generators.marginal_cost
    dynamic = n.get_switchable_as_dense("Generator", "marginal_cost")

    marginal_cost = dynamic + co2_cost.reindex(columns=dynamic.columns, fill_value=0)
    n.generators_t.marginal_cost = marginal_cost.loc[:, marginal_cost.ne(static).any()]


def set_line_s_max_pu(n, s_max_pu=0.7):
    n.lines["s_max_pu"] = s_max_pu
    logger.info(f"N-1 security margin of lines set to {s_max_pu}")


def set_transmission_limit(n, kind, factor, costs, Nyears=1):
    links_dc_b = n.links.carrier == "DC" if not n.links.empty else pd.Series()

    _lines_s_nom = (
        np.sqrt(3)
        * n.lines.type.map(n.line_types.i_nom)
        * n.lines.num_parallel
        * n.lines.bus0.map(n.buses.v_nom)
    )
    lines_s_nom = n.lines.s_nom.where(n.lines.type == "", _lines_s_nom)

    col = "capital_cost" if kind == "c" else "length"
    ref = (
        lines_s_nom @ n.lines[col]
        + n.links.loc[links_dc_b, "p_nom"] @ n.links.loc[links_dc_b, col]
    )

    set_transmission_costs(n, costs)

    if factor == "opt" or float(factor) > 1.0:
        n.lines["s_nom_min"] = lines_s_nom
        n.lines["s_nom_extendable"] = True

        n.links.loc[links_dc_b, "p_nom_min"] = n.links.loc[links_dc_b, "p_nom"]
        n.links.loc[links_dc_b, "p_nom_extendable"] = True

    if factor != "opt":
        con_type = "expansion_cost" if kind == "c" else "volume_expansion"
        rhs = float(factor) * ref
        n.add(
            "GlobalConstraint",
            f"l{kind}_limit",
            type=f"transmission_{con_type}_limit",
            sense="<=",
            constant=rhs,
            carrier_attribute="AC, DC",
        )

    return n


def average_every_nhours(n, offset, drop_leap_day=False):
    logger.info(f"Resampling the network to {offset}")
    m = n.copy(snapshots=[])

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    sns = snapshot_weightings.index
    if drop_leap_day:
        sns = sns[~((sns.month == 2) & (sns.day == 29))]
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                pnl[k] = df.resample(offset).mean()

    return m


def apply_time_segmentation(n, segments, solver_name="cbc"):
    logger.info(f"Aggregating time series to {segments} segments.")
    try:
        import tsam.timeseriesaggregation as tsam
    except ImportError:
        raise ModuleNotFoundError(
            "Optional dependency 'tsam' not found.Install via 'pip install tsam'"
        )

    p_max_pu_norm = n.generators_t.p_max_pu.max()
    p_max_pu = n.generators_t.p_max_pu / p_max_pu_norm

    load_norm = n.loads_t.p_set.max()
    load = n.loads_t.p_set / load_norm

    inflow_norm = n.storage_units_t.inflow.max()
    inflow = n.storage_units_t.inflow / inflow_norm

    raw = pd.concat([p_max_pu, load, inflow], axis=1, sort=False)

    agg = tsam.TimeSeriesAggregation(
        raw,
        hoursPerPeriod=len(raw),
        noTypicalPeriods=1,
        noSegments=int(segments),
        segmentation=True,
        solver=solver_name,
    )

    segmented = agg.createTypicalPeriods()

    weightings = segmented.index.get_level_values("Segment Duration")
    offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
    snapshots = [n.snapshots[0] + pd.Timedelta(f"{offset}h") for offset in offsets]

    n.set_snapshots(pd.DatetimeIndex(snapshots, name="name"))
    n.snapshot_weightings = pd.Series(
        weightings, index=snapshots, name="weightings", dtype="float64"
    )

    segmented.index = snapshots
    n.generators_t.p_max_pu = segmented[n.generators_t.p_max_pu.columns] * p_max_pu_norm
    n.loads_t.p_set = segmented[n.loads_t.p_set.columns] * load_norm
    n.storage_units_t.inflow = segmented[n.storage_units_t.inflow.columns] * inflow_norm

    return n


def enforce_autarky(n, only_crossborder=False):
    if only_crossborder:
        lines_rm = n.lines.loc[
            n.lines.bus0.map(n.buses.country) != n.lines.bus1.map(n.buses.country)
        ].index
        links_rm = n.links.loc[
            n.links.bus0.map(n.buses.country) != n.links.bus1.map(n.buses.country)
        ].index
    else:
        lines_rm = n.lines.index
        links_rm = n.links.loc[n.links.carrier == "DC"].index
    n.remove("Line", lines_rm)
    n.remove("Link", links_rm)


def cap_transmission_capacity(
    n,
    line_max=None,
    link_max=None,
    line_max_extension=None,
    link_max_extension=None,
    line_max_pu=None,
    link_max_pu=None,
):
    """
    Cap transmission capacity for AC lines and DC links.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    line_max : float, optional
        Absolute upper limit for AC line capacity [MW]. If None, no limit is applied.
    link_max : float, optional
        Absolute upper limit for DC link capacity [MW]. If None, no limit is applied.
    line_max_extension : float, optional
        Maximum extension per AC line [MW]. If None, no limit is applied.
    link_max_extension : float, optional
        Maximum extension per DC link [MW]. If None, no limit is applied.
    line_max_pu : float, optional
        Set N-1 security margin for AC lines (e.g., 0.7 for 70% utilization).
        If None, s_max_pu is not modified.
    link_max_pu : float, optional
        Set maximum utilization for DC links (e.g., 0.7 for 70% utilization).
        If None, p_max_pu is not modified.

    Notes
    -----
    All parameters accept None to skip that particular constraint. This allows
    selective application of limits without needing to specify all parameters.
    """
    # Set N-1 security margin (s_max_pu) for AC lines if specified
    if line_max_pu is not None:
        n.lines["s_max_pu"] = line_max_pu
        logger.info(f"N-1 security margin of lines set to {line_max_pu}")

    # Set maximum utilization (p_max_pu) for DC links if specified
    if link_max_pu is not None:
        hvdc = n.links.index[n.links.carrier == "DC"]
        n.links.loc[hvdc, "p_max_pu"] = link_max_pu
        logger.info(f"Maximum utilization of DC links set to {link_max_pu}")

    # Apply line capacity extension limit if specified
    if (
        line_max_extension is not None
        and np.isfinite(line_max_extension)
        and line_max_extension > 0
    ):
        logger.info(f"Limiting AC line extensions to {line_max_extension} MW")
        n.lines["s_nom_max"] = n.lines["s_nom"] + line_max_extension

    # Apply link capacity extension limit if specified
    if (
        link_max_extension is not None
        and np.isfinite(link_max_extension)
        and link_max_extension > 0
    ):
        logger.info(f"Limiting DC link extensions to {link_max_extension} MW")
        hvdc = n.links.index[n.links.carrier == "DC"]
        n.links.loc[hvdc, "p_nom_max"] = n.links.loc[hvdc, "p_nom"] + link_max_extension

    # Apply absolute line capacity limit if specified
    if line_max is not None and np.isfinite(line_max):
        n.lines["s_nom_max"] = n.lines.s_nom_max.clip(upper=line_max)

    # Apply absolute link capacity limit if specified
    if link_max is not None and np.isfinite(link_max):
        n.links["p_nom_max"] = n.links.p_nom_max.clip(upper=link_max)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_network",
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    logger.warning(
        "Deprecated: final network preparation is handled by "
        "compose_network.py under 'NETWORK PREPARATION (from prepare_network.py)' "
        "(see prepare_network_for_solving). Use compose_network.py instead."
    )

    n = pypsa.Network(snakemake.input[0])
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0
    costs = load_costs(snakemake.input.costs)

    set_line_s_max_pu(n, snakemake.params.lines["s_max_pu"])

    # temporal averaging
    time_resolution = snakemake.params.time_resolution
    is_string = isinstance(time_resolution, str)
    if is_string and time_resolution.lower().endswith("h"):
        n = average_every_nhours(n, time_resolution, snakemake.params.drop_leap_day)

    # segments with package tsam
    if is_string and time_resolution.lower().endswith("seg"):
        solver_name = snakemake.config["solving"]["solver"]["name"]
        segments = int(time_resolution.replace("seg", ""))
        n = apply_time_segmentation(n, segments, solver_name)

    if snakemake.params.co2limit_enable:
        add_co2limit(n, snakemake.params.co2limit, None, Nyears)

    if snakemake.params.gaslimit_enable:
        add_gaslimit(n, snakemake.params.gaslimit, Nyears)

    maybe_adjust_costs_and_potentials(n, snakemake.params["adjustments"])

    emission_prices = snakemake.params.emission_prices
    if emission_prices["co2_monthly_prices"]:
        logger.info(
            "Setting time dependent emission prices according spot market price"
        )
        add_dynamic_emission_prices(n, snakemake.input.co2_price)
    elif emission_prices["enable"]:
        if isinstance(emission_prices["co2"], dict):
            logger.warning(
                "Not setting emission prices on generators and storage units, "
                "due to their configuration per planning horizon"
            )
        elif isinstance(emission_prices["co2"], float):
            add_emission_prices(n, dict(co2=emission_prices["co2"]))

    kind = snakemake.params.transmission_limit[0]
    factor = snakemake.params.transmission_limit[1:]
    set_transmission_limit(n, kind, factor, costs, Nyears)

    cap_transmission_capacity(
        n,
        line_max=snakemake.params.lines.get("s_nom_max", np.inf),
        link_max=snakemake.params.links.get("p_nom_max", np.inf),
        line_max_extension=snakemake.params.lines.get("max_extension", np.inf),
        link_max_extension=snakemake.params.links.get("max_extension", np.inf),
    )

    if snakemake.params.autarky["enable"]:
        only_crossborder = snakemake.params.autarky["by_country"]
        enforce_autarky(n, only_crossborder=only_crossborder)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
