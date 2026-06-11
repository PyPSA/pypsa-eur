# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves optimal operation and capacity for a network with the option to
iteratively optimize while updating line reactances.

This script is used for optimizing the electrical network as well as the
sector coupled network.

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.

The optimization is based on the :func:`network.optimize` function.
Additionally, some extra constraints specified in :mod:`solve_network` are added.

.. note::

    The rules ``solve_elec_networks`` and ``solve_sector_networks`` run
    the workflow for all scenarios in the configuration file (``scenario:``)
    based on the rule :mod:`solve_network`.
"""

import importlib
import logging
import os
import re
import sys
from functools import partial
from pathlib import Path
from typing import Any

import linopy
import linopy.io
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
import yaml
from linopy.remote.oetc import OetcCredentials, OetcHandler, OetcSettings
from pypsa.descriptors import get_activity_mask
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

from scripts._benchmark import memory_logger
from scripts._helpers import (
    PYPSA_V1,
    configure_logging,
    get,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_brownfield import disable_grid_expansion_if_limit_hit

logger = logging.getLogger(__name__)


DEFAULT_DIAGNOSTIC_SLACK_PATTERNS = {
    "global": [r"^GlobalConstraint-", r"co2_sequestration_limit"],
    "cdr": [r"^cdr_", r"^CDR-", r"^CO2AnnualCDRSeq"],
    "biomass": [r"biomass", r"biogas", r"solid.?biomass"],
    "stores": [r"^Store-.*(energy_balance|e_cyclic|e_initial|e_lower|e_upper)"],
    "load_balance": [r"^Bus-nodal_balance$"],
    "imports": [r"^import_limit$"],
    "growth": [r"growth", r"agg_p_nom", r"bau_mincaps", r"safe_mintotalcap"],
}


def _transmission_expansion_cost_reference(n: pypsa.Network) -> float:
    """Return the AC/DC reference cost used for cost-based grid limits."""
    links_dc_b = n.links.carrier == "DC" if not n.links.empty else pd.Series()

    lines_s_nom = n.lines.s_nom
    typed_lines = n.lines.type != ""
    if typed_lines.any():
        typed_s_nom = (
            np.sqrt(3)
            * n.lines.loc[typed_lines, "type"].map(n.line_types.i_nom)
            * n.lines.loc[typed_lines, "num_parallel"]
            * n.lines.loc[typed_lines, "bus0"].map(n.buses.v_nom)
        )
        lines_s_nom = lines_s_nom.where(~typed_lines, typed_s_nom)

    return (
        lines_s_nom @ n.lines.capital_cost
        + n.links.loc[links_dc_b, "p_nom"] @ n.links.loc[links_dc_b, "capital_cost"]
    )


def _rescale_transmission_expansion_cost_limits(
    n: pypsa.Network, reference_before: float
) -> None:
    """Keep cost-based grid limits consistent after line/link cost perturbations."""
    if reference_before <= 0:
        return

    glcs = n.global_constraints.query("type == 'transmission_expansion_cost_limit'")
    if glcs.empty:
        return

    reference_after = _transmission_expansion_cost_reference(n)
    if reference_after <= 0:
        return

    factor = reference_after / reference_before
    for name in glcs.index:
        n.global_constraints.at[name, "constant"] *= factor
        logger.info(
            "Rescaled transmission expansion cost limit %s after noisy_costs "
            "(factor %.9f)",
            name,
            factor,
        )


def patch_linopy_highspy_name_export() -> None:
    """Avoid allocating millions of string labels when exporting MPS via HiGHS."""
    if getattr(linopy.io, "_pypsa_highspy_name_export_patched", False):
        return

    def _to_highspy_without_names(m, explicit_coordinate_names: bool = False):
        if m.variables.sos:
            raise NotImplementedError(
                "SOS constraints are not supported by the HiGHS direct API. "
                "Use io_api='lp' instead."
            )

        import highspy
        from scipy.sparse import triu

        M = m.matrices
        h = highspy.Highs()
        h.addVars(len(M.vlabels), M.lb, M.ub)
        if len(m.binaries) + len(m.integers):
            vtypes = M.vtypes
            labels = np.arange(len(vtypes))[(vtypes == "B") | (vtypes == "I")]
            n = len(labels)
            h.changeColsIntegrality(n, labels, np.ones_like(labels))
            if len(m.binaries):
                labels = np.arange(len(vtypes))[vtypes == "B"]
                n = len(labels)
                h.changeColsBounds(
                    n, labels, np.zeros_like(labels), np.ones_like(labels)
                )

        h.changeColsCost(len(M.c), np.arange(len(M.c), dtype=np.int32), M.c)

        A = M.A
        if A is not None:
            A = A.tocsr()
            num_cons = A.shape[0]
            lower = np.where(M.sense != "<", M.b, -np.inf)
            upper = np.where(M.sense != ">", M.b, np.inf)
            h.addRows(num_cons, lower, upper, A.nnz, A.indptr, A.indices, A.data)

        # Skip row/column names entirely to avoid large string allocations.
        h.passModel(h.getLp())

        Q = M.Q
        if Q is not None:
            Q = triu(Q)
            Q = Q.tocsr()
            num_vars = Q.shape[0]
            h.passHessian(num_vars, Q.nnz, 1, Q.indptr, Q.indices, Q.data)

        if m.objective.sense == "max":
            h.changeObjectiveSense(highspy.ObjSense.kMaximize)

        return h

    linopy.io.to_highspy = _to_highspy_without_names
    linopy.model.to_highspy = _to_highspy_without_names
    linopy.model.Model.to_highspy = _to_highspy_without_names
    linopy.io._pypsa_highspy_name_export_patched = True


# Allow for PyPSA versions <0.35
if PYPSA_V1:
    pypsa.network.power_flow.logger.setLevel(logging.WARNING)
else:
    pypsa.pf.logger.setLevel(logging.WARNING)


class ObjectiveValueError(Exception):
    pass


def add_land_use_constraint_perfect(n: pypsa.Network) -> None:
    """
    Add global constraints for tech capacity limit.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance

    Returns
    -------
    pypsa.Network
        Network with added land use constraints
    """
    logger.info("Add land-use constraint for perfect foresight")

    def compress_series(s):
        def process_group(group):
            if group.nunique() == 1:
                return pd.Series(group.iloc[0], index=[None])
            else:
                return group

        return s.groupby(level=[0, 1]).apply(process_group)

    def new_index_name(t):
        # Convert all elements to string and filter out None values
        parts = [str(x) for x in t if x is not None]
        # Join with space, but use a dash for the last item if not None
        return " ".join(parts[:2]) + (f"-{parts[-1]}" if len(parts) > 2 else "")

    def check_p_min_p_max(p_nom_max):
        p_nom_min = n.generators[ext_i].groupby(grouper).sum().p_nom_min
        p_nom_min = p_nom_min.reindex(p_nom_max.index)
        check = (
            p_nom_min.groupby(level=[0, 1]).sum()
            > p_nom_max.groupby(level=[0, 1]).min()
        )
        if check.sum():
            logger.warning(
                f"summed p_min_pu values at node larger than technical potential {check[check].index}"
            )

    grouper = [n.generators.carrier, n.generators.bus, n.generators.build_year]
    ext_i = n.generators.p_nom_extendable
    # get technical limit per node and investment period
    p_nom_max = n.generators[ext_i].groupby(grouper).min().p_nom_max
    # drop carriers without tech limit
    p_nom_max = p_nom_max[~p_nom_max.isin([np.inf, np.nan])]
    # carrier
    carriers = p_nom_max.index.get_level_values(0).unique()
    gen_i = n.generators[(n.generators.carrier.isin(carriers)) & (ext_i)].index
    n.generators.loc[gen_i, "p_nom_min"] = 0
    # check minimum capacities
    check_p_min_p_max(p_nom_max)
    # drop multi entries in case p_nom_max stays constant in different periods
    # p_nom_max = compress_series(p_nom_max)
    # adjust name to fit syntax of nominal constraint per bus
    df = p_nom_max.reset_index()
    df["name"] = df.apply(
        lambda row: f"nom_max_{row['carrier']}"
        + (f"_{row['build_year']}" if row["build_year"] is not None else ""),
        axis=1,
    )

    for name in df.name.unique():
        df_carrier = df[df.name == name]
        bus = df_carrier.bus
        n.buses.loc[bus, name] = df_carrier.p_nom_max.values


def add_land_use_constraint(n: pypsa.Network, planning_horizons: str) -> None:
    """
    Add land use constraints for renewable energy potential.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    planning_horizons : str
        The planning horizon year as string

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with constraints added
    """
    # warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in [
        "solar",
        "solar rooftop",
        "solar-hsat",
        "onwind",
        "offwind-ac",
        "offwind-dc",
        "offwind-float",
    ]:
        ext_i = (n.generators.carrier == carrier) & ~n.generators.p_nom_extendable
        grouper = n.generators.loc[ext_i].index.str.replace(
            f" {carrier}.*$", "", regex=True
        )
        existing = n.generators.loc[ext_i, "p_nom"].groupby(grouper).sum()
        existing.index += f" {carrier}-{planning_horizons}"
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    # check if existing capacities are larger than technical potential
    existing_large = n.generators[
        n.generators["p_nom_min"] > n.generators["p_nom_max"]
    ].index
    if len(existing_large):
        logger.warning(
            f"Existing capacities larger than technical potential for {existing_large},\
                        adjust technical potential to existing capacities"
        )
        n.generators.loc[existing_large, "p_nom_max"] = n.generators.loc[
            existing_large, "p_nom_min"
        ]

    n.generators["p_nom_max"] = n.generators["p_nom_max"].clip(lower=0)


def add_solar_potential_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add constraint to make sure the sum capacity of all solar technologies (fixed, tracking, ets. ) is below the region potential.

    Example:
    ES1 0: total solar potential is 10 GW, meaning:
           solar potential : 10 GW
           solar-hsat potential : 8 GW (solar with single axis tracking is assumed to have higher land use)
    The constraint ensures that:
           solar_p_nom + solar_hsat_p_nom * 1.13 <= 10 GW
    """
    land_use_factors = {
        "solar-hsat": config["renewable"]["solar"]["capacity_per_sqkm"]
        / config["renewable"]["solar-hsat"]["capacity_per_sqkm"],
    }
    rename = {} if PYPSA_V1 else {"Generator-ext": "Generator"}

    solar_carriers = ["solar", "solar-hsat"]
    solar = n.generators[
        n.generators.carrier.isin(solar_carriers) & n.generators.p_nom_extendable
    ].index

    solar_today = n.generators[
        (n.generators.carrier == "solar") & (n.generators.p_nom_extendable)
    ].index
    solar_hsat = n.generators[(n.generators.carrier == "solar-hsat")].index

    if solar.empty:
        return

    land_use = pd.DataFrame(1, index=solar, columns=["land_use_factor"])
    for carrier, factor in land_use_factors.items():
        land_use = land_use.apply(
            lambda x: (x * factor) if carrier in x.name else x, axis=1
        )

    location = pd.Series(n.buses.index, index=n.buses.index)
    ggrouper = n.generators.loc[solar].bus
    rhs = (
        n.generators.loc[solar_today, "p_nom_max"]
        .groupby(n.generators.loc[solar_today].bus.map(location))
        .sum()
        - n.generators.loc[solar_hsat, "p_nom"]
        .groupby(n.generators.loc[solar_hsat].bus.map(location))
        .sum()
        * land_use_factors["solar-hsat"]
    ).clip(lower=0)

    lhs = (
        (n.model["Generator-p_nom"].rename(rename).loc[solar] * land_use.squeeze())
        .groupby(ggrouper)
        .sum()
    )

    logger.info("Adding solar potential constraint.")
    n.model.add_constraints(lhs <= rhs, name="solar_potential")


def add_co2_sequestration_limit(
    n: pypsa.Network,
    limit_dict: dict[str, float],
    planning_horizons: str | None,
) -> None:
    """
    Add a global constraint on the amount of Mt CO2 that can be sequestered.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    limit_dict : dict[str, float]
        CO2 sequestration potential limit constraints by year.
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    """

    if not n.investment_periods.empty:
        nyears = n.snapshot_weightings.groupby(level="period").generators.sum() / 8760
        periods = n.investment_periods
        limit = pd.Series(
            {period: nyears[period] * get(limit_dict, period) for period in periods}
        )
        limit.index = limit.index.map(lambda s: f"co2_sequestration_limit-{s}")
        names = limit.index
    else:
        nyears = n.snapshot_weightings.generators.sum() / 8760
        limit = get(limit_dict, int(planning_horizons)) * nyears
        periods = np.nan
        names = "co2_sequestration_limit"

    n.add(
        "GlobalConstraint",
        names,
        sense=">=",
        constant=-limit * 1e6,
        type="operational_limit",
        carrier_attribute="co2 sequestered",
        investment_period=periods,
    )


CDR_CREDIT_ORIGINS = ("dac", "biogenic", "fossil")
CDR_INFRASTRUCTURE_PREFIXES = (
    "co2 pipeline",
    "co2 sequestered",
    "co2 vent",
    "co2 provenance",
)
BIOGENIC_TOKENS = ("biomass", "biogas", "biosng", "btl", "fuelwood", "msw", "bio")


def _eligible_cdr_scopes(cdr_credit_scope: str | list[str]) -> set[str]:
    if cdr_credit_scope == "all_sequestration":
        return set(CDR_CREDIT_ORIGINS)
    if isinstance(cdr_credit_scope, str):
        scopes = {cdr_credit_scope}
    else:
        scopes = set(cdr_credit_scope)
    scopes.discard("all_sequestration")
    return scopes


def _classify_co2_origin(carrier_name: str) -> str:
    carrier = str(carrier_name).lower()
    if "dac" in carrier:
        return "dac"
    if any(token in carrier for token in BIOGENIC_TOKENS):
        return "biogenic"
    return "fossil"


def _get_cdr_credit_prices_for_period(
    cdr_credit_price,
    cdr_credit_scope: str | list[str],
    cdr_credit_prices_by_scope,
    planning_horizons: str | None,
) -> dict[str, float]:
    period = int(planning_horizons) if planning_horizons is not None else None
    if cdr_credit_prices_by_scope:
        return {
            str(scope): float(get(traj, period))
            for scope, traj in cdr_credit_prices_by_scope.items()
            if float(get(traj, period)) != 0.0
        }

    base_price = float(get(cdr_credit_price, period))
    if base_price == 0.0:
        return {}

    scopes = _eligible_cdr_scopes(cdr_credit_scope)
    if not scopes:
        return {}
    return {scope: base_price for scope in scopes}


def _co2_stored_buses(n: pypsa.Network) -> pd.Index:
    return pd.Index(
        n.buses.index[n.buses.carrier.astype(str).str.startswith("co2 stored")]
    )


def _capture_term_data(n: pypsa.Network) -> pd.DataFrame:
    co2_stored_buses = set(_co2_stored_buses(n))
    rows = []
    for link_name in n.links.index:
        carrier = str(n.links.at[link_name, "carrier"]).lower()
        if carrier.startswith(CDR_INFRASTRUCTURE_PREFIXES):
            continue

        origin = _classify_co2_origin(n.links.at[link_name, "carrier"])
        for idx in [1, 2, 3, 4]:
            bus_col = f"bus{idx}"
            eff_col = f"efficiency{idx}" if idx > 1 else "efficiency"
            if bus_col in n.links.columns and eff_col in n.links.columns:
                bus = n.links.at[link_name, bus_col]
                eff = n.links.at[link_name, eff_col]
                if pd.notna(bus) and pd.notna(eff):
                    eff = float(eff)
                    if bus in co2_stored_buses and eff > 0:
                        rows.append(
                            {
                                "link": link_name,
                                "bus": bus,
                                "origin": origin,
                                "coeff": eff,
                            }
                        )

    if not rows:
        return pd.DataFrame(columns=["link", "bus", "origin", "coeff"])
    return pd.DataFrame(rows)


def _capture_link_data(n: pypsa.Network) -> pd.DataFrame:
    co2_stored_buses = set(_co2_stored_buses(n))
    rows = []
    for link_name in n.links.index:
        carrier = str(n.links.at[link_name, "carrier"]).lower()
        if carrier.startswith(CDR_INFRASTRUCTURE_PREFIXES):
            continue
        coeff = 0.0
        atm_withdrawal = 0.0
        for idx in [1, 2, 3, 4]:
            bus_col = f"bus{idx}"
            eff_col = f"efficiency{idx}" if idx > 1 else "efficiency"
            if bus_col in n.links.columns and eff_col in n.links.columns:
                bus = n.links.at[link_name, bus_col]
                eff = n.links.at[link_name, eff_col]
                if pd.notna(bus) and pd.notna(eff):
                    if bus in co2_stored_buses and float(eff) > 0:
                        coeff += float(eff)
                    elif bus == "co2 atmosphere" and float(eff) < 0:
                        atm_withdrawal += -float(eff)
        if coeff > 0.0:
            rows.append(
                {
                    "link": link_name,
                    "origin": _classify_co2_origin(n.links.at[link_name, "carrier"]),
                    "coeff": coeff,
                    "atm_withdrawal": atm_withdrawal,
                }
            )
    return pd.DataFrame(rows)


def _withdrawal_term_data(n: pypsa.Network) -> pd.DataFrame:
    co2_stored_buses = set(_co2_stored_buses(n))
    rows = []
    for link_name in n.links.index:
        carrier = str(n.links.at[link_name, "carrier"]).lower()
        if carrier.startswith(CDR_INFRASTRUCTURE_PREFIXES):
            continue

        for idx in [1, 2, 3, 4]:
            bus_col = f"bus{idx}"
            eff_col = f"efficiency{idx}" if idx > 1 else "efficiency"
            if bus_col in n.links.columns and eff_col in n.links.columns:
                bus = n.links.at[link_name, bus_col]
                eff = n.links.at[link_name, eff_col]
                if pd.notna(bus) and pd.notna(eff):
                    eff = float(eff)
                    if bus in co2_stored_buses and eff < 0:
                        rows.append(
                            {
                                "term": f"co2_use_term_{len(rows)}",
                                "link": link_name,
                                "bus": bus,
                                "coeff": -eff,
                            }
                        )

    if not rows:
        return pd.DataFrame(columns=["term", "link", "bus", "coeff"])
    return pd.DataFrame(rows).set_index("term", drop=False)


def _co2_pipeline_links(n: pypsa.Network) -> pd.Index:
    return n.links.index[
        n.links.carrier.astype(str).str.lower().str.startswith("co2 pipeline")
    ]


def _sequestration_links(n: pypsa.Network) -> pd.Index:
    return n.links.index[n.links.carrier == "co2 sequestered"]


def _vent_links(n: pypsa.Network) -> pd.Index:
    return n.links.index[n.links.carrier == "co2 vent"]


def _previous_snapshots(snapshot_index: pd.Index) -> pd.Index:
    if snapshot_index.empty:
        return snapshot_index

    prev_positions = np.empty(len(snapshot_index), dtype=int)
    if isinstance(snapshot_index, pd.MultiIndex):
        if "period" in snapshot_index.names:
            period_values = snapshot_index.get_level_values("period")
        else:
            period_values = snapshot_index.get_level_values(0)

        for period in pd.Index(period_values).unique():
            positions = np.flatnonzero(period_values == period)
            prev_positions[positions] = np.roll(positions, 1)
    else:
        prev_positions[:] = np.roll(np.arange(len(snapshot_index)), 1)

    return snapshot_index[prev_positions]


def _link_term_timeseries(
    link_p: xr.DataArray,
    link_dim: str,
    terms: pd.DataFrame,
) -> xr.DataArray | float:
    if terms.empty:
        return 0.0

    coeffs = terms.groupby("link")["coeff"].sum()
    dispatch = link_p.sel({link_dim: coeffs.index.tolist()})
    coeffs_da = xr.DataArray(
        coeffs.to_numpy(),
        dims=[link_dim],
        coords={link_dim: coeffs.index.tolist()},
    )
    return (dispatch * coeffs_da).sum(link_dim)


def _weighted_link_term_total(
    link_p: xr.DataArray,
    link_dim: str,
    snap_dim: str,
    snapshot_index: pd.Index,
    snapshot_weightings: pd.Series,
    terms: pd.DataFrame,
) -> xr.DataArray | float:
    if terms.empty:
        return 0.0

    coeffs = terms.groupby("link")["coeff"].sum()
    dispatch = link_p.sel({link_dim: coeffs.index.tolist()})
    coeffs_da = xr.DataArray(
        np.outer(snapshot_weightings.to_numpy(), coeffs.to_numpy()),
        dims=[snap_dim, link_dim],
        coords={snap_dim: snapshot_index, link_dim: coeffs.index.tolist()},
    )
    return (dispatch * coeffs_da).sum()


def _sanitize_constraint_key(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(value)).strip("_")


def add_cdr_credit_accounting(
    n: pypsa.Network,
    planning_horizons: str | None,
) -> None:
    cdr_credit_timing = n.params.get("cdr_credit_timing", "capture")
    if cdr_credit_timing != "sequestration" or planning_horizons is None:
        return

    cdr_credit_scope = n.params.get("cdr_credit_scope") or []
    eligible_scopes = _eligible_cdr_scopes(cdr_credit_scope)
    cdr_credit_standalone = bool(n.params.get("cdr_credit_standalone", False))
    standalone_origins = {"dac", "biogenic"} if cdr_credit_standalone else set()
    accounting_origins = eligible_scopes | standalone_origins

    prices = _get_cdr_credit_prices_for_period(
        cdr_credit_price=n.params.get("cdr_credit_price", 0.0),
        cdr_credit_scope=cdr_credit_scope,
        cdr_credit_prices_by_scope=n.params.get("cdr_credit_prices_by_scope", {}),
        planning_horizons=planning_horizons,
    )
    cdr_credit_limit = n.params.get("cdr_credit_limit")
    cdr_credit_limit_by_year = n.params.get("cdr_credit_limit_by_year")

    # Provide explicit finite upper bounds so Gurobi's barrier doesn't treat
    # these variables as free.  Unbounded CDR variables poison the homogeneous
    # barrier on 2040+ brownfield models (infeasible_or_unbounded with code 4).
    _limit_src = cdr_credit_limit_by_year or cdr_credit_limit
    if _limit_src and planning_horizons:
        nyears = n.snapshot_weightings.generators.sum() / 8760
        _cdr_var_ub = get(_limit_src, int(planning_horizons)) * nyears * 1e6
    else:
        _cdr_var_ub = np.inf

    needs_accounting = (
        bool(prices)
        or bool(cdr_credit_limit)
        or bool(cdr_credit_limit_by_year)
        or cdr_credit_standalone
    )
    if not needs_accounting or not accounting_origins:
        return

    capture = _capture_link_data(n)
    relevant_capture = capture[capture.origin.isin(accounting_origins)]
    eligible_capture = capture[capture.origin.isin(eligible_scopes)]
    if relevant_capture.empty:
        logger.info("add_cdr_credit_accounting: no relevant capture links found, skipping.")
        return

    capture_terms = _capture_term_data(n)
    if capture_terms.empty:
        logger.info("add_cdr_credit_accounting: no CO2 capture terms found, skipping.")
        return

    sequest_links = _sequestration_links(n)
    if sequest_links.empty:
        logger.info("add_cdr_credit_accounting: no physical sequestration links found, skipping.")
        return

    link_p = n.model["Link-p"]
    snap_dim = link_p.dims[0]
    link_dim = link_p.dims[1]
    snapshot_index = link_p.coords[snap_dim].to_index()
    generator_weightings = n.snapshot_weightings.generators.reindex(snapshot_index)
    store_weightings = n.snapshot_weightings.stores.reindex(snapshot_index)

    if generator_weightings.isna().any() or store_weightings.isna().any():
        raise ValueError("Snapshot weightings could not be aligned with solver snapshots.")

    credited_index = pd.Index(sorted(eligible_scopes), name="cdr_origin")
    credited = None
    if not credited_index.empty:
        n.model.add_variables(0, _cdr_var_ub, coords=[credited_index], name="CDR-credited")
        credited = n.model["CDR-credited"]

    # -----------------------------------------------------------------------
    # ANNUAL CDR SEQUESTRATION TRACKING
    # Replaces per-snapshot CO2 origin inventory tracking (~1 M variables)
    # with annual-level attribution variables (one scalar per accounting
    # origin, typically 2: dac + biogenic).
    #
    # The per-snapshot approach caused Gurobi numerical trouble on 2040
    # brownfield models (dual bound exceeds primal after ~170 barrier
    # iterations).  Annual tracking preserves all economically relevant
    # constraints:
    #   - CDR credits <= annual sequestration attributed to eligible origin
    #   - CDR credits <= annual capture by eligible origin
    #   - total CDR attribution <= total physical annual sequestration
    #   - ETS cancellation applied to DAC's annual attributed sequestration
    # -----------------------------------------------------------------------

    generator_weightings_da = xr.DataArray(
        generator_weightings.to_numpy(),
        dims=[snap_dim],
        coords={snap_dim: snapshot_index},
    )

    # Annual total physical sequestration across all links (tCO2/yr)
    annual_total_seq = (
        link_p.sel({link_dim: sequest_links.tolist()}) * generator_weightings_da
    ).sum()

    accounting_index = pd.Index(sorted(accounting_origins), name="cdr_origin")
    n.model.add_variables(
        0,
        _cdr_var_ub,
        coords=[accounting_index],
        name="CO2AnnualCDRSeq",
    )
    annual_cdr_seq = n.model["CO2AnnualCDRSeq"]

    # Total CDR attribution <= total physical annual sequestration
    n.model.add_constraints(
        annual_cdr_seq.sum("cdr_origin") <= annual_total_seq,
        name="cdr_annual_seq_total_limit",
    )

    # Per-origin attribution <= annual capture by that origin
    for origin in sorted(accounting_origins):
        origin_capture_terms = capture_terms[capture_terms.origin == origin]
        if origin_capture_terms.empty:
            n.model.add_constraints(
                annual_cdr_seq.sel({"cdr_origin": origin}) <= 0.0,
                name=f"cdr_annual_seq_capture_limit-{_sanitize_constraint_key(origin)}",
            )
        else:
            annual_capture = _weighted_link_term_total(
                link_p=link_p,
                link_dim=link_dim,
                snap_dim=snap_dim,
                snapshot_index=snapshot_index,
                snapshot_weightings=generator_weightings,
                terms=origin_capture_terms,
            )
            n.model.add_constraints(
                annual_cdr_seq.sel({"cdr_origin": origin}) <= annual_capture,
                name=f"cdr_annual_seq_capture_limit-{_sanitize_constraint_key(origin)}",
            )

    # sequestered_by_origin: linopy variable per origin, used for CDR credit
    # limits and ETS cancellation in the objective below.
    sequestered_by_origin = {
        origin: annual_cdr_seq.sel({"cdr_origin": origin})
        for origin in accounting_origins
    }

    if credited is not None:
        for origin in credited_index:
            origin_capture = eligible_capture[eligible_capture.origin == origin]
            if origin_capture.empty:
                n.model.add_constraints(
                    credited.loc[origin] <= 0.0,
                    name=f"cdr_credited_capture_limit-{origin}",
                )
                continue
            captured_expr = _weighted_link_term_total(
                link_p=link_p,
                link_dim=link_dim,
                snap_dim=snap_dim,
                snapshot_index=snapshot_index,
                snapshot_weightings=generator_weightings,
                terms=capture_terms[capture_terms.origin == origin],
            )
            n.model.add_constraints(
                credited.loc[origin] <= captured_expr,
                name=f"cdr_credited_capture_limit-{origin}",
            )
            n.model.add_constraints(
                credited.loc[origin] <= sequestered_by_origin[origin],
                name=f"cdr_credited_sequestration_limit-{origin}",
            )

    # Build objective adjustment: CDR credit revenue minus ETS cancellation for
    # sequestered CO2 only.  Under standalone mode, DAC/BECCS should not earn
    # both ETS and CDR credit.  By deferring the ETS cancellation to here (rather
    # than as a fixed marginal cost at build time) we cancel ETS only for the
    # CO2 that actually reaches geological sequestration.  CO2 diverted to CCU
    # retains its ETS credit at capture and pays ETS at combustion, keeping
    # that loop ETS-neutral.
    objective_delta = None

    if prices and credited is not None:
        revenue = sum(
            float(prices.get(origin, 0.0)) * credited.loc[origin]
            for origin in credited_index
        )
        objective_delta = -revenue

    # ETS cancellation for standalone DAC/BECCS is applied unconditionally at
    # build time in prepare_sector_network.py (apply_cdr_credit_to_eligible_capture_links).
    # Doing it there avoids the optimizer driving CO2AnnualCDRSeq to zero to dodge
    # the net cost (ets_price - credit_price) when ets_price > credit_price.

    if objective_delta is not None:
        objective = n.model.objective.expression + objective_delta
        n.model.add_objective(objective, overwrite=True, sense=n.model.objective.sense)


def add_cdr_credit_limit(
    n: pypsa.Network,
    limit_dict: dict[str, float],
    planning_horizons: str | None,
    cdr_credit_scope: str | list[str],
    cdr_credit_timing: str = "capture",
) -> None:
    """Add a constraint capping the total CDR credits issued per year (Mt CO2)."""
    eligible_scopes = _eligible_cdr_scopes(cdr_credit_scope)
    if not eligible_scopes:
        return

    nyears = n.snapshot_weightings.generators.sum() / 8760
    limit = get(limit_dict, int(planning_horizons)) * nyears

    if cdr_credit_timing == "sequestration":
        if "CDR-credited" not in n.model.variables:
            logger.info("add_cdr_credit_limit: credited CDR variables missing, skipping.")
            return
        credited = n.model["CDR-credited"]
        lhs = credited.sum()
        n.model.add_constraints(lhs <= limit * 1e6, name="cdr_credit_limit")
        return

    dim = "name" if PYPSA_V1 else "Link"
    link_p = n.model["Link-p"]
    weightings = n.snapshot_weightings.generators
    snap_dim = link_p.dims[0]
    snap_coords = link_p.coords[snap_dim].values

    capture = _capture_link_data(n)
    capture = capture[capture.origin.isin(eligible_scopes)]
    if capture.empty:
        logger.info("add_cdr_credit_limit: no eligible CDR links found, skipping.")
        return

    eligible_links = capture.link.tolist()
    coeffs_da = xr.DataArray(
        np.outer(weightings.values, capture.coeff.to_numpy()),
        dims=[snap_dim, dim],
        coords={snap_dim: snap_coords, dim: eligible_links},
    )
    lhs = (link_p.sel({dim: eligible_links}) * coeffs_da).sum()
    n.model.add_constraints(lhs <= limit * 1e6, name="cdr_credit_limit")


def _series_from_solution(
    variable: linopy.variables.Variable | None,
    model: linopy.Model | None = None,
    name: str | None = None,
) -> pd.Series:
    solutions = []
    if variable is not None:
        for attr in ("solution", "sol"):
            try:
                solution = getattr(variable, attr, None)
            except Exception as exc:
                logger.debug("Could not read %s for variable %s: %s", attr, name, exc)
                solution = None
            if solution is not None:
                solutions.append(solution)

    if model is not None and name is not None:
        try:
            solution = getattr(model, "solution", None)
            if solution is not None and name in solution:
                solutions.append(solution[name])
        except Exception as exc:
            logger.debug("Could not read %s from model.solution: %s", name, exc)

    for solution in solutions:
        if hasattr(solution, "to_series"):
            series = solution.to_series()
        else:
            series = pd.Series(solution)

        if isinstance(series.index, pd.MultiIndex):
            if len(series.index.names) == 1:
                series.index = series.index.get_level_values(0)

        series = pd.to_numeric(series, errors="coerce").dropna()
        if not series.empty:
            return series.astype(float)

    return pd.Series(dtype=float)


def _model_variable(n: pypsa.Network, name: str) -> linopy.variables.Variable | None:
    if not hasattr(n, "model"):
        return None
    try:
        return n.model[name]
    except Exception:
        return None


def _annual_capture_proxy_by_origin(n: pypsa.Network) -> dict[str, float]:
    capture = _capture_link_data(n)
    if capture.empty:
        return {origin: 0.0 for origin in CDR_CREDIT_ORIGINS}

    weights = n.snapshot_weightings.generators
    out = {origin: 0.0 for origin in CDR_CREDIT_ORIGINS}

    for origin in out:
        origin_capture = capture[capture.origin == origin]
        if origin_capture.empty:
            continue
        coeffs = origin_capture.groupby("link")["coeff"].sum()
        dispatch = n.links_t.p0[coeffs.index]
        weighted = dispatch.mul(weights, axis=0).mul(coeffs, axis=1)
        out[origin] = float(weighted.to_numpy().sum())

    return out


def _annual_physical_sequestration(n: pypsa.Network) -> float:
    links = _sequestration_links(n)
    if links.empty:
        return 0.0

    weights = n.snapshot_weightings.generators
    dispatch = n.links_t.p0[links]
    return float(dispatch.mul(weights, axis=0).to_numpy().sum())


def export_cdr_credit_accounting(
    n: pypsa.Network,
    output_path: str | Path,
    planning_horizons: str | None,
) -> None:
    period = int(planning_horizons) if planning_horizons is not None else None
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    prices = _get_cdr_credit_prices_for_period(
        cdr_credit_price=n.params.get("cdr_credit_price", 0.0),
        cdr_credit_scope=n.params.get("cdr_credit_scope") or [],
        cdr_credit_prices_by_scope=n.params.get("cdr_credit_prices_by_scope", {}),
        planning_horizons=planning_horizons,
    )
    eligible_scopes = sorted(_eligible_cdr_scopes(n.params.get("cdr_credit_scope") or []))
    limit_dict = n.params.get("cdr_credit_limit_by_year") or n.params.get("cdr_credit_limit")
    credit_limit_mt = float(get(limit_dict, period)) if (limit_dict and period is not None) else np.nan

    capture_proxy = _annual_capture_proxy_by_origin(n)
    physical_seq_t = _annual_physical_sequestration(n)

    model = getattr(n, "model", None)
    credited_series = _series_from_solution(
        _model_variable(n, "CDR-credited"),
        model=model,
        name="CDR-credited",
    )
    attributed_series = _series_from_solution(
        _model_variable(n, "CO2AnnualCDRSeq"),
        model=model,
        name="CO2AnnualCDRSeq",
    )

    use_solver_export = not credited_series.empty
    cdr_credit_timing = n.params.get("cdr_credit_timing", "capture")
    needs_sequestration_accounting = (
        cdr_credit_timing == "sequestration"
        and bool(eligible_scopes)
        and (
            bool(prices)
            or bool(n.params.get("cdr_credit_limit"))
            or bool(n.params.get("cdr_credit_limit_by_year"))
            or bool(n.params.get("cdr_credit_standalone", False))
        )
    )

    solver_values_missing = needs_sequestration_accounting and (
        credited_series.empty or attributed_series.empty
    )
    if solver_values_missing:
        logger.warning(
            "CDR accounting is configured for sequestration timing, but solver "
            "values for CDR-credited/CO2AnnualCDRSeq are unavailable. "
            "Falling back to capture_proxy for credited amounts."
        )

    def _mt(series: pd.Series, key: str) -> float:
        if key not in series.index:
            return 0.0
        return float(series.loc[key]) / 1e6

    if solver_values_missing:
        # ── Deterministic sequestration proration ────────────────────────────
        # The "co2 stored" bus is a single fungible commodity: DAC, biogenic and
        # fossil CO2 are commingled before geological sequestration, and synfuels
        # (Fischer-Tropsch / Sabatier / methanolisation) draw back out of the same
        # pool.  There is therefore no physical fact about which captured tonne is
        # the one that ends up underground, which is exactly why the solver-side
        # CO2AnnualCDRSeq attribution variable is underdetermined (the optimiser
        # leaves it at an arbitrary feasible point / its solution is not reliably
        # recoverable post-solve -> attributed_series is empty here every time).
        #
        # With cdr_credit_timing == "sequestration" we instead allocate the *actual*
        # physical geological sequestration to capture origins pro-rata by captured
        # volume, and credit only the eligible (DAC + biogenic) shares.  Fossil CO2
        # is stored but earns no CDR credit.  By construction this guarantees
        #   credited_origin   <= capture_origin           (seq <= total capture), and
        #   credited_total    <= physical_sequestration,
        # so the headline credited CDR can never exceed what is physically stored.
        cap_dac = capture_proxy.get("dac", 0.0)
        cap_biogenic = capture_proxy.get("biogenic", 0.0)
        cap_fossil = capture_proxy.get("fossil", 0.0)
        cap_total = cap_dac + cap_biogenic + cap_fossil
        if cap_total > 0 and physical_seq_t > 0:
            credited_dac_mt = (physical_seq_t * cap_dac / cap_total) / 1e6
            credited_biogenic_mt = (physical_seq_t * cap_biogenic / cap_total) / 1e6
            # Respect an explicit (binding) credit limit if one is configured.
            credited_total_tmp = credited_dac_mt + credited_biogenic_mt
            if (
                not np.isnan(credit_limit_mt)
                and credited_total_tmp > credit_limit_mt
                and credited_total_tmp > 0
            ):
                scale = credit_limit_mt / credited_total_tmp
                credited_dac_mt *= scale
                credited_biogenic_mt *= scale
        else:
            # Nothing eligible captured, or nothing physically sequestered.
            credited_dac_mt = 0.0
            credited_biogenic_mt = 0.0
        method = "sequestration_proration"
        solver_values_missing = False  # credited values are now available
        logger.info(
            "CDR sequestration proration: credited dac=%.2f Mt, biogenic=%.2f Mt "
            "(physical_seq=%.2f Mt; capture dac/bio/fossil=%.2f/%.2f/%.2f Mt)",
            credited_dac_mt,
            credited_biogenic_mt,
            physical_seq_t / 1e6,
            cap_dac / 1e6,
            cap_biogenic / 1e6,
            cap_fossil / 1e6,
        )
    elif use_solver_export:
        credited_dac_mt = _mt(credited_series, "dac")
        credited_biogenic_mt = _mt(credited_series, "biogenic")
        method = "solver_export"
    else:
        credited_dac_mt = capture_proxy["dac"] / 1e6
        credited_biogenic_mt = capture_proxy["biogenic"] / 1e6
        method = "capture_proxy"

    row = {
        "planning_horizon": period,
        "cdr_credit_timing": cdr_credit_timing,
        "cdr_credit_standalone": bool(n.params.get("cdr_credit_standalone", False)),
        "method": method,
        "accounting_error": (
            "missing_solver_values_for_CDR-credited_or_CO2AnnualCDRSeq"
            if solver_values_missing
            else ""
        ),
        "eligible_scopes": ",".join(eligible_scopes),
        "price_dac_eur_per_tco2": float(prices.get("dac", 0.0)),
        "price_biogenic_eur_per_tco2": float(prices.get("biogenic", 0.0)),
        "credit_limit_mtco2_per_yr": credit_limit_mt,
        "credited_dac_mtco2_per_yr": credited_dac_mt,
        "credited_biogenic_mtco2_per_yr": credited_biogenic_mt,
        "credited_total_mtco2_per_yr": credited_dac_mt + credited_biogenic_mt,
        "attributed_dac_mtco2_per_yr": _mt(attributed_series, "dac"),
        "attributed_biogenic_mtco2_per_yr": _mt(attributed_series, "biogenic"),
        "attributed_total_mtco2_per_yr": (
            (_mt(attributed_series, "dac") + _mt(attributed_series, "biogenic"))
            if not attributed_series.empty
            else np.nan
        ),
        "capture_proxy_dac_mtco2_per_yr": capture_proxy["dac"] / 1e6,
        "capture_proxy_biogenic_mtco2_per_yr": capture_proxy["biogenic"] / 1e6,
        "capture_proxy_total_mtco2_per_yr": (
            capture_proxy["dac"] + capture_proxy["biogenic"]
        )
        / 1e6,
        "physical_sequestration_mtco2_per_yr": physical_seq_t / 1e6,
    }
    tolerance_mt = 1e-3
    credited_total = row["credited_total_mtco2_per_yr"]
    attributed_total = row["attributed_total_mtco2_per_yr"]
    physical_seq_mt = row["physical_sequestration_mtco2_per_yr"]
    capture_proxy_total = row["capture_proxy_total_mtco2_per_yr"]

    row["valid_credited_within_limit"] = False if solver_values_missing else (
        True
        if np.isnan(credit_limit_mt)
        else credited_total <= credit_limit_mt + tolerance_mt
    )
    row["valid_credited_within_capture_proxy"] = (
        False
        if solver_values_missing
        else credited_total <= capture_proxy_total + tolerance_mt
    )
    row["valid_credited_within_attributed"] = (
        False
        if np.isnan(attributed_total)
        else credited_total <= attributed_total + tolerance_mt
    )
    # Primary invariant for the sequestration_proration method: the headline
    # credited CDR must never exceed what is physically sequestered.
    row["valid_credited_within_physical_sequestration"] = (
        False
        if solver_values_missing
        else credited_total <= physical_seq_mt + tolerance_mt
    )
    row["valid_attributed_within_physical_sequestration"] = (
        False
        if np.isnan(attributed_total)
        else attributed_total <= physical_seq_mt + tolerance_mt
    )

    pd.DataFrame([row]).to_csv(output_path, index=False)
    logger.info("Exported CDR credit accounting to %s", output_path.resolve())


def export_cdr_credit_accounting_failure(
    n: pypsa.Network,
    output_path: str | Path,
    planning_horizons: str | None,
    error: Exception,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "planning_horizon": (
                    int(planning_horizons) if planning_horizons is not None else np.nan
                ),
                "cdr_credit_timing": n.params.get("cdr_credit_timing", "capture"),
                "cdr_credit_standalone": bool(
                    n.params.get("cdr_credit_standalone", False)
                ),
                "method": "export_failed",
                "accounting_error": str(error),
                "valid_credited_within_limit": False,
                "valid_credited_within_capture_proxy": False,
                "valid_credited_within_attributed": False,
                "valid_attributed_within_physical_sequestration": False,
            }
        ]
    ).to_csv(output_path, index=False)
    logger.warning(
        "Exported flagged CDR credit accounting failure row to %s",
        output_path.resolve(),
    )


def add_carbon_constraint(n: pypsa.Network, snapshots: pd.DatetimeIndex) -> None:
    glcs = n.global_constraints.query('type == "co2_atmosphere"')
    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last = n.snapshot_weightings.reset_index().groupby("period").last()
            last_i = last.set_index([last.index, last.timestep]).index
            final_e = n.model["Store-e"].loc[last_i, stores.index]
            time_valid = int(glc.loc["investment_period"])
            time_i = pd.IndexSlice[time_valid, :]
            lhs = final_e.loc[time_i, :] - final_e.shift(snapshot=1).loc[time_i, :]

            rhs = glc.constant
            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def add_carbon_budget_constraint(n: pypsa.Network, snapshots: pd.DatetimeIndex) -> None:
    glcs = n.global_constraints.query('type == "Co2Budget"')
    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last = n.snapshot_weightings.reset_index().groupby("period").last()
            last_i = last.set_index([last.index, last.timestep]).index
            final_e = n.model["Store-e"].loc[last_i, stores.index]
            time_valid = int(glc.loc["investment_period"])
            time_i = pd.IndexSlice[time_valid, :]
            weighting = n.investment_period_weightings.loc[time_valid, "years"]
            lhs = final_e.loc[time_i, :] * weighting

            rhs = glc.constant
            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def add_max_growth(n: pypsa.Network, opts: dict) -> None:
    """
    Add maximum growth rates for different carriers.
    """

    # take maximum yearly difference between investment periods since historic growth is per year
    factor = n.investment_period_weightings.years.max() * opts["factor"]
    for carrier in opts["max_growth"].keys():
        max_per_period = opts["max_growth"][carrier] * factor
        logger.info(
            f"set maximum growth rate per investment period of {carrier} to {max_per_period} GW."
        )
        n.carriers.loc[carrier, "max_growth"] = max_per_period * 1e3

    for carrier in opts["max_relative_growth"].keys():
        max_r_per_period = opts["max_relative_growth"][carrier]
        logger.info(
            f"set maximum relative growth per investment period of {carrier} to {max_r_per_period}."
        )
        n.carriers.loc[carrier, "max_relative_growth"] = max_r_per_period


def _myopic_period_years(n: pypsa.Network, planning_horizons: str) -> int:
    horizons = sorted(
        int(year)
        for year in n.config.get("scenario", {}).get("planning_horizons", [])
    )
    period = int(planning_horizons)

    if period in horizons:
        index = horizons.index(period)
        if index > 0:
            return period - horizons[index - 1]
        if len(horizons) > 1:
            return horizons[1] - horizons[0]

    return 10


def add_max_growth_myopic(
    n: pypsa.Network, opts: dict, planning_horizons: str
) -> None:
    """
    Apply per-carrier maximum growth constraints for myopic optimisation.

    The perfect-foresight version (add_max_growth) sets n.carriers.max_growth
    which PyPSA enforces automatically via investment_period_weightings.  In
    myopic mode investment_period_weightings is empty, so those carrier
    attributes are never applied.  This function adds equivalent explicit
    linopy constraints:

        sum(p_nom_opt) <= existing_mw + rate_GW_yr * period_years * factor

    and, when max_relative_growth is configured:

        sum(extendable p_nom_opt) <= max_relative * brownfield_mw

    where brownfield_mw is the non-extendable (already-built) fleet inherited
    from the previous horizon.  The relative constraint is skipped when the
    brownfield fleet is zero to avoid forcing new-technology deployment to zero.
    """
    if not opts.get("enable", False):
        return

    period_years = _myopic_period_years(n, planning_horizons)
    factor = period_years * opts["factor"]

    for carrier, rate_gw_yr in opts["max_growth"].items():
        max_new_mw = rate_gw_yr * factor * 1e3  # GW/yr to MW

        gen_i = n.generators.index[
            (n.generators.carrier == carrier) & n.generators.p_nom_extendable
        ]
        link_i = n.links.index[
            (n.links.carrier == carrier) & n.links.p_nom_extendable
        ]

        if gen_i.empty and link_i.empty:
            logger.debug(
                "max_growth_myopic: no extendable components for %s, skipping.",
                carrier,
            )
            continue

        existing_mw = n.generators.loc[gen_i, "p_nom"].sum() + n.links.loc[link_i, "p_nom"].sum()
        cap_mw = existing_mw + max_new_mw

        lhs_parts = []
        if not gen_i.empty:
            lhs_parts.append(n.model["Generator-p_nom"].loc[gen_i].sum())
        if not link_i.empty:
            lhs_parts.append(n.model["Link-p_nom"].loc[link_i].sum())
        lhs = lhs_parts[0] if len(lhs_parts) == 1 else sum(lhs_parts[1:], lhs_parts[0])

        key = _sanitize_constraint_key(carrier)
        n.model.add_constraints(lhs <= cap_mw, name=f"max_growth_myopic_{key}")
        logger.info(
            f"max_growth_myopic: {carrier} total cap = {cap_mw/1e3:.0f} GW "
            f"(existing {existing_mw/1e3:.0f} GW + new <= {max_new_mw/1e3:.0f} GW "
            f"= {rate_gw_yr} GW/yr x {period_years} yr x {opts['factor']})"
        )

    for carrier, max_relative in opts.get("max_relative_growth", {}).items():
        gen_ext_i = n.generators.index[
            (n.generators.carrier == carrier) & n.generators.p_nom_extendable
        ]
        link_ext_i = n.links.index[
            (n.links.carrier == carrier) & n.links.p_nom_extendable
        ]

        if gen_ext_i.empty and link_ext_i.empty:
            logger.debug(
                "max_growth_myopic (relative): no extendable components for %s, skipping.",
                carrier,
            )
            continue

        gen_fixed_i = n.generators.index[
            (n.generators.carrier == carrier) & ~n.generators.p_nom_extendable
        ]
        link_fixed_i = n.links.index[
            (n.links.carrier == carrier) & ~n.links.p_nom_extendable
        ]
        brownfield_mw = n.generators.loc[gen_fixed_i, "p_nom"].sum() + n.links.loc[
            link_fixed_i, "p_nom"
        ].sum()

        if brownfield_mw == 0:
            logger.debug(
                "max_growth_myopic (relative): brownfield fleet for %s is zero, skipping "
                "relative cap to avoid forcing new technology to zero.",
                carrier,
            )
            continue

        cap_mw = max_relative * brownfield_mw

        lhs_parts = []
        if not gen_ext_i.empty:
            lhs_parts.append(n.model["Generator-p_nom"].loc[gen_ext_i].sum())
        if not link_ext_i.empty:
            lhs_parts.append(n.model["Link-p_nom"].loc[link_ext_i].sum())
        lhs = lhs_parts[0] if len(lhs_parts) == 1 else sum(lhs_parts[1:], lhs_parts[0])

        key = _sanitize_constraint_key(carrier)
        n.model.add_constraints(lhs <= cap_mw, name=f"max_growth_myopic_{key}_relative")
        logger.info(
            f"max_growth_myopic (relative): {carrier} new cap <= {cap_mw/1e3:.0f} GW "
            f"({max_relative}x brownfield {brownfield_mw/1e3:.0f} GW)"
        )


def add_retrofit_gas_boiler_constraint(
    n: pypsa.Network, snapshots: pd.DatetimeIndex
) -> None:
    """
    Allow retrofitting of existing gas boilers to H2 boilers and impose load-following must-run condition on existing gas boilers.
    Modifies the network in place, no return value.

    n : pypsa.Network
        The PyPSA network to be modified
    snapshots : pd.DatetimeIndex
        The snapshots of the network
    """
    c = "Link"
    logger.info("Add constraint for retrofitting gas boilers to H2 boilers.")
    # existing gas boilers
    mask = n.links.carrier.str.contains("gas boiler") & ~n.links.p_nom_extendable
    gas_i = n.links[mask].index
    mask = n.links.carrier.str.contains("retrofitted H2 boiler")
    h2_i = n.links[mask].index

    n.links.loc[gas_i, "p_nom_extendable"] = True
    p_nom = n.links.loc[gas_i, "p_nom"]
    n.links.loc[gas_i, "p_nom"] = 0

    # heat profile
    cols = n.loads_t.p_set.columns[
        n.loads_t.p_set.columns.str.contains("heat")
        & ~n.loads_t.p_set.columns.str.contains("industry")
        & ~n.loads_t.p_set.columns.str.contains("agriculture")
    ]
    profile = n.loads_t.p_set[cols].div(
        n.loads_t.p_set[cols].groupby(level=0).max(), level=0
    )
    # to deal if max value is zero
    profile.fillna(0, inplace=True)
    profile.rename(columns=n.loads.bus.to_dict(), inplace=True)
    profile = profile.reindex(columns=n.links.loc[gas_i, "bus1"])
    profile.columns = gas_i

    rhs = profile.mul(p_nom)

    dispatch = n.model["Link-p"]
    active = get_activity_mask(n, c, snapshots, gas_i)
    rhs = rhs[active]
    if PYPSA_V1:
        p_gas = dispatch.sel(name=gas_i)
        p_h2 = dispatch.sel(name=h2_i)
    else:
        p_gas = dispatch.sel(Link=gas_i)
        p_h2 = dispatch.sel(Link=h2_i)

    lhs = p_gas + p_h2

    n.model.add_constraints(lhs == rhs, name="gas_retrofit")


def prepare_network(
    n: pypsa.Network,
    solve_opts: dict,
    foresight: str,
    planning_horizons: str | None,
    co2_sequestration_potential: dict[str, float],
    limit_max_growth: dict[str, Any] | None = None,
) -> None:
    """
    Prepare network with various constraints and modifications.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    solve_opts : Dict
        Dictionary of solving options containing clip_p_max_pu, load_shedding etc.
    foresight : str
        Planning foresight type ('myopic' or 'perfect')
    planning_horizons : str or None
        The current planning horizon year or None for perfect foresight
    co2_sequestration_potential : Dict[str, float]
        CO2 sequestration potential constraints by year

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with added constraints
    """
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.links_t.p_max_pu,
            n.links_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df.abs() > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    if load_shedding := solve_opts.get("load_shedding"):
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        n.add("Carrier", "load")
        buses_i = n.buses.index
        if isinstance(load_shedding, bool):
            load_shedding = 1e5  # Eur/MWh

        n.add(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            marginal_cost=load_shedding,  # Eur/MWh
            p_nom=np.inf,
        )

    if solve_opts.get("curtailment_mode"):
        n.add("Carrier", "curtailment", color="#fedfed", nice_name="Curtailment")
        n.generators_t.p_min_pu = n.generators_t.p_max_pu
        buses_i = n.buses.query("carrier == 'AC'").index
        n.add(
            "Generator",
            buses_i,
            suffix=" curtailment",
            bus=buses_i,
            p_min_pu=-1,
            p_max_pu=0,
            marginal_cost=-0.1,
            carrier="curtailment",
            p_nom=1e6,
        )

    if solve_opts.get("noisy_costs"):
        transmission_reference_before = _transmission_expansion_cost_reference(n)

        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

        _rescale_transmission_expansion_cost_limits(
            n, transmission_reference_before
        )
        disable_grid_expansion_if_limit_hit(n)

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if foresight == "myopic" and planning_horizons:
        add_land_use_constraint(n, planning_horizons)

    if foresight == "perfect":
        add_land_use_constraint_perfect(n)
        if limit_max_growth is not None and limit_max_growth["enable"]:
            add_max_growth(n, limit_max_growth)

    if n.stores.carrier.eq("co2 sequestered").any():
        limit_dict = co2_sequestration_potential
        add_co2_sequestration_limit(
            n, limit_dict=limit_dict, planning_horizons=planning_horizons
        )

def add_CCL_constraints(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> None:
    """
    Add CCL (country & carrier limit) constraint to the network.

    Add minimum and maximum levels of generator nominal capacity per carrier
    for individual countries. Opts and path for agg_p_nom_minmax.csv must be defined
    in config.yaml. Default file is available at data/agg_p_nom_minmax.csv.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : dict
        Configuration dictionary
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Example
    -------
    scenario:
        opts: [Co2L-CCL-24h]
    electricity:
        agg_p_nom_limits: data/agg_p_nom_minmax.csv
    """

    assert planning_horizons is not None, (
        "add_CCL_constraints are not implemented for perfect foresight, yet"
    )

    agg_p_nom_minmax = pd.read_csv(
        config["solving"]["agg_p_nom_limits"]["file"], index_col=[0, 1], header=[0, 1]
    )[planning_horizons]
    logger.info("Adding generation capacity constraints per carrier and country")
    p_nom = n.model["Generator-p_nom"]

    gens = n.generators.query("p_nom_extendable")

    if not PYPSA_V1:
        gens = gens.rename_axis(index="Generator-ext")

    if config["solving"]["agg_p_nom_limits"]["agg_offwind"]:
        rename_offwind = {
            "offwind-ac": "offwind-all",
            "offwind-dc": "offwind-all",
            "offwind-float": "offwind-all",
            "offwind": "offwind-all",
        }
        gens = gens.replace(rename_offwind)
    if config["solving"]["agg_p_nom_limits"]["agg_solar"]:
        rename_solar = {
            "solar": "solar-all",
            "solar-hsat": "solar-all",
            "solar rooftop": "solar-all",
        }
        gens = gens.replace(rename_solar)
    grouper = pd.concat([gens.bus.map(n.buses.country), gens.carrier], axis=1)
    lhs = p_nom.groupby(grouper).sum().rename(bus="country")

    if config["solving"]["agg_p_nom_limits"]["include_existing"]:
        gens_cst = n.generators.query("~p_nom_extendable").rename_axis(
            index="Generator-cst"
        )
        gens_cst = gens_cst[
            (gens_cst["build_year"] + gens_cst["lifetime"]) >= int(planning_horizons)
        ]
        if config["solving"]["agg_p_nom_limits"]["agg_offwind"]:
            gens_cst = gens_cst.replace(rename_offwind)
        if config["solving"]["agg_p_nom_limits"]["agg_solar"]:
            gens_cst = gens_cst.replace(rename_solar)
        rhs_cst = (
            pd.concat(
                [gens_cst.bus.map(n.buses.country), gens_cst[["carrier", "p_nom"]]],
                axis=1,
            )
            .groupby(["bus", "carrier"])
            .sum()
        )
        rhs_cst.index = rhs_cst.index.rename({"bus": "country"})
        rhs_min = agg_p_nom_minmax["min"].dropna()
        idx_min = rhs_min.index.join(rhs_cst.index, how="left")
        rhs_min = rhs_min.reindex(idx_min).fillna(0)
        rhs = (rhs_min - rhs_cst.reindex(idx_min).fillna(0).p_nom).dropna()
        rhs[rhs < 0] = 0
        minimum = xr.DataArray(rhs).rename(dim_0="group")
    else:
        minimum = xr.DataArray(agg_p_nom_minmax["min"].dropna()).rename(dim_0="group")

    index = minimum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) >= minimum.loc[index], name="agg_p_nom_min"
        )

    if config["solving"]["agg_p_nom_limits"]["include_existing"]:
        rhs_max = agg_p_nom_minmax["max"].dropna()
        idx_max = rhs_max.index.join(rhs_cst.index, how="left")
        rhs_max = rhs_max.reindex(idx_max).fillna(0)
        rhs = (rhs_max - rhs_cst.reindex(idx_max).fillna(0).p_nom).dropna()
        rhs[rhs < 0] = 0
        maximum = xr.DataArray(rhs).rename(dim_0="group")
    else:
        maximum = xr.DataArray(agg_p_nom_minmax["max"].dropna()).rename(dim_0="group")

    index = maximum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) <= maximum.loc[index], name="agg_p_nom_max"
        )


def add_EQ_constraints(n, o, scaling=1e-1):
    """
    Add equity constraints to the network.

    Currently this is only implemented for the electricity sector only.

    Opts must be specified in the config.yaml.

    Parameters
    ----------
    n : pypsa.Network
    o : str

    Example
    -------
    scenario:
        opts: [Co2L-EQ0.7-24h]

    Require each country or node to on average produce a minimal share
    of its total electricity consumption itself. Example: EQ0.7c demands each country
    to produce on average at least 70% of its consumption; EQ0.7 demands
    each node to produce on average at least 70% of its consumption.
    """
    # TODO: Generalize to cover myopic and other sectors?
    float_regex = r"[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])
    if o[-1] == "c":
        ggrouper = n.generators.bus.map(n.buses.country)
        lgrouper = n.loads.bus.map(n.buses.country)
        sgrouper = n.storage_units.bus.map(n.buses.country)
    else:
        ggrouper = n.generators.bus
        lgrouper = n.loads.bus
        sgrouper = n.storage_units.bus
    load = (
        n.snapshot_weightings.generators
        @ n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    )
    inflow = (
        n.snapshot_weightings.stores
        @ n.storage_units_t.inflow.groupby(sgrouper, axis=1).sum()
    )
    inflow = inflow.reindex(load.index).fillna(0.0)
    rhs = scaling * (level * load - inflow)
    p = n.model["Generator-p"]
    lhs_gen = (
        (p * (n.snapshot_weightings.generators * scaling))
        .groupby(ggrouper.to_xarray())
        .sum()
        .sum("snapshot")
    )
    # TODO: double check that this is really needed, why do have to subtract the spillage
    if not n.storage_units_t.inflow.empty:
        spillage = n.model["StorageUnit-spill"]
        lhs_spill = (
            (spillage * (-n.snapshot_weightings.stores * scaling))
            .groupby(sgrouper.to_xarray())
            .sum()
            .sum("snapshot")
        )
        lhs = lhs_gen + lhs_spill
    else:
        lhs = lhs_gen
    n.model.add_constraints(lhs >= rhs, name="equity_min")


def add_BAU_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add business-as-usual (BAU) constraints for minimum capacities.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network instance
    config : dict
        Configuration dictionary containing BAU minimum capacities
    """
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    p_nom = n.model["Generator-p_nom"]
    ext_i = n.generators.query("p_nom_extendable")
    ext_carrier_i = xr.DataArray(ext_i.carrier)
    if not PYPSA_V1:
        ext_carrier_i = ext_carrier_i.rename_axis("Generator-ext")
    lhs = p_nom.groupby(ext_carrier_i).sum()
    rhs = mincaps[lhs.indexes["carrier"]].rename_axis("carrier")
    n.model.add_constraints(lhs >= rhs, name="bau_mincaps")


# TODO: think about removing or make per country
def add_SAFE_constraints(n, config):
    """
    Add a capacity reserve margin of a certain fraction above the peak demand.
    Renewable generators and storage do not contribute. Ignores network.

    Parameters
    ----------
        n : pypsa.Network
        config : dict

    Example
    -------
    config.yaml requires to specify opts:

    scenario:
        opts: [Co2L-SAFE-24h]
    electricity:
        SAFE_reservemargin: 0.1
    Which sets a reserve margin of 10% above the peak demand.
    """
    peakdemand = n.loads_t.p_set.sum(axis=1).max()
    margin = 1.0 + config["electricity"]["SAFE_reservemargin"]
    reserve_margin = peakdemand * margin
    conventional_carriers = config["electricity"]["conventional_carriers"]  # noqa: F841
    ext_gens_i = n.generators.query(
        "carrier in @conventional_carriers & p_nom_extendable"
    ).index
    p_nom = n.model["Generator-p_nom"].loc[ext_gens_i]
    lhs = p_nom.sum()
    exist_conv_caps = n.generators.query(
        "~p_nom_extendable & carrier in @conventional_carriers"
    ).p_nom.sum()
    rhs = reserve_margin - exist_conv_caps
    n.model.add_constraints(lhs >= rhs, name="safe_mintotalcap")


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.

    Parameters
    ----------
        n : pypsa.Network
        sns: pd.DatetimeIndex
        config : dict

    Example:
    --------
    config.yaml requires to specify operational_reserve:
    operational_reserve: # like https://genxproject.github.io/GenX/dev/core/#Reserves
        activate: true
        epsilon_load: 0.02 # percentage of load at each snapshot
        epsilon_vres: 0.02 # percentage of VRES at each snapshot
        contingency: 400000 # MW
    """
    reserve_config = config["electricity"]["operational_reserve"]
    EPSILON_LOAD = reserve_config["epsilon_load"]
    EPSILON_VRES = reserve_config["epsilon_vres"]
    CONTINGENCY = reserve_config["contingency"]

    # Reserve Variables
    n.model.add_variables(
        0, np.inf, coords=[sns, n.generators.index], name="Generator-r"
    )
    reserve = n.model["Generator-r"]
    summed_reserve = reserve.sum("Generator")

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        p_nom_vres = n.model["Generator-p_nom"].loc[vres_i.intersection(ext_i)]
        if not PYPSA_V1:
            p_nom_vres = p_nom_vres.rename({"Generator-ext": "Generator"})
        lhs = summed_reserve + (
            p_nom_vres * (-EPSILON_VRES * xr.DataArray(capacity_factor))
        ).sum("Generator")

        # Total demand per t
        demand = get_as_dense(n, "Load", "p_set").sum(axis=1)

        # VRES potential of non extendable generators
        capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
        renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
        potential = (capacity_factor * renewable_capacity).sum(axis=1)

        # Right-hand-side
        rhs = EPSILON_LOAD * demand + EPSILON_VRES * potential + CONTINGENCY

        n.model.add_constraints(lhs >= rhs, name="reserve_margin")

    # additional constraint that capacity is not exceeded
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = n.model["Generator-p"]
    reserve = n.model["Generator-r"]

    capacity_variable = n.model["Generator-p_nom"]
    if not PYPSA_V1:
        capacity_variable = capacity_variable.rename({"Generator-ext": "Generator"})
    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    lhs = dispatch + reserve - capacity_variable * xr.DataArray(p_max_pu[ext_i])

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)

    n.model.add_constraints(lhs <= rhs, name="Generator-p-reserve-upper")


def add_TES_energy_to_power_ratio_constraints(n: pypsa.Network) -> None:
    """
    Add TES constraints to the network.

    For each TES storage unit, enforce:
        Store-e_nom - etpr * Link-p_nom == 0

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA network with TES and heating sectors enabled.

    Raises
    ------
    ValueError
        If no valid TES storage or charger links are found.
    RuntimeError
        If the TES storage and charger indices do not align.
    """
    indices_charger_p_nom_extendable = n.links.index[
        n.links.index.str.contains("water tanks charger|water pits charger")
        & n.links.p_nom_extendable
    ]
    indices_stores_e_nom_extendable = n.stores.index[
        n.stores.index.str.contains("water tanks|water pits")
        & n.stores.e_nom_extendable
    ]

    if indices_charger_p_nom_extendable.empty or indices_stores_e_nom_extendable.empty:
        logger.warning(
            "No valid extendable charger links or stores found for TES energy-to-power constraints.Not enforcing TES energy-to-power ratio constraints!"
        )
        return

    energy_to_power_ratio_values = n.links.loc[
        indices_charger_p_nom_extendable, "energy to power ratio"
    ].values

    linear_expr_list = []
    for charger, tes, energy_to_power_value in zip(
        indices_charger_p_nom_extendable,
        indices_stores_e_nom_extendable,
        energy_to_power_ratio_values,
    ):
        charger_var = n.model["Link-p_nom"].loc[charger]
        if not tes == charger.replace(" charger", ""):
            # e.g. "DE0 0 urban central water tanks charger-2050" -> "DE0 0 urban central water tanks-2050"
            raise RuntimeError(
                f"Charger {charger} and TES {tes} do not match. "
                "Ensure that the charger and TES are in the same location and refer to the same technology."
            )
        store_var = n.model["Store-e_nom"].loc[tes]
        linear_expr = store_var - energy_to_power_value * charger_var
        linear_expr_list.append(linear_expr)

    # Merge the individual expressions
    dim = "Store-ext, Link-ext" if PYPSA_V1 else "name"
    merged_expr = linopy.expressions.merge(
        linear_expr_list, dim=dim, cls=type(linear_expr_list[0])
    )

    n.model.add_constraints(merged_expr == 0, name="TES_energy_to_power_ratio")


def add_TES_charger_ratio_constraints(n: pypsa.Network) -> None:
    """
    Add TES charger ratio constraints.

    For each TES unit, enforce:
        Link-p_nom(charger) - efficiency * Link-p_nom(discharger) == 0

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA network with TES and heating sectors enabled.

    Raises
    ------
    ValueError
        If no valid TES discharger or charger links are found.
    RuntimeError
        If the charger and discharger indices do not align.
    """
    indices_charger_p_nom_extendable = n.links.index[
        n.links.index.str.contains(
            "water tanks charger|water pits charger|aquifer thermal energy storage charger"
        )
        & n.links.p_nom_extendable
    ]
    indices_discharger_p_nom_extendable = n.links.index[
        n.links.index.str.contains(
            "water tanks discharger|water pits discharger|aquifer thermal energy storage discharger"
        )
        & n.links.p_nom_extendable
    ]

    if (
        indices_charger_p_nom_extendable.empty
        or indices_discharger_p_nom_extendable.empty
    ):
        logger.warning(
            "No valid extendable TES discharger or charger links found for TES charger ratio constraints. Not enforcing TES charger_ratio constraints."
        )
        return

    for charger, discharger in zip(
        indices_charger_p_nom_extendable, indices_discharger_p_nom_extendable
    ):
        if not charger.replace(" charger", " ") == discharger.replace(
            " discharger", " "
        ):
            # e.g. "DE0 0 urban central water tanks charger-2050" -> "DE0 0 urban central water tanks-2050"
            raise RuntimeError(
                f"Charger {charger} and discharger {discharger} do not match. "
                "Ensure that the charger and discharger are in the same location and refer to the same technology."
            )

    eff_discharger = n.links.efficiency[indices_discharger_p_nom_extendable].values
    lhs = (
        n.model["Link-p_nom"].loc[indices_charger_p_nom_extendable]
        - n.model["Link-p_nom"].loc[indices_discharger_p_nom_extendable]
        * eff_discharger
    )

    n.model.add_constraints(lhs == 0, name="TES_charger_ratio")


def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.index.str.contains("battery discharger")
    charger_bool = n.links.index.str.contains("battery charger")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


def add_lossy_bidirectional_link_constraints(n):
    if not n.links.p_nom_extendable.any() or not any(n.links.get("reversed", [])):
        return

    carriers = n.links.loc[n.links.reversed, "carrier"].unique()  # noqa: F841
    backwards = n.links.query(
        "carrier in @carriers and p_nom_extendable and reversed"
    ).index
    forwards = backwards.str.replace("-reversed", "")
    lhs = n.model["Link-p_nom"].loc[backwards]
    rhs = n.model["Link-p_nom"].loc[forwards]
    n.model.add_constraints(lhs == rhs, name="Link-bidirectional_sync")


def add_chp_constraints(n):
    electric = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("electric")
    )
    heat = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("heat")
    )

    electric_ext = n.links[electric].query("p_nom_extendable").index
    heat_ext = n.links[heat].query("p_nom_extendable").index

    electric_fix = n.links[electric].query("~p_nom_extendable").index
    heat_fix = n.links[heat].query("~p_nom_extendable").index

    p = n.model["Link-p"]  # dimension: [time, link]

    # output ratio between heat and electricity and top_iso_fuel_line for extendable
    if not electric_ext.empty:
        p_nom = n.model["Link-p_nom"]

        lhs = (
            p_nom.loc[electric_ext]
            * (n.links.p_nom_ratio * n.links.efficiency)[electric_ext].values
            - p_nom.loc[heat_ext] * n.links.efficiency[heat_ext].values
        )
        n.model.add_constraints(lhs == 0, name="chplink-fix_p_nom_ratio")

        rename = {} if PYPSA_V1 else {"Link-ext": "Link"}
        lhs = (
            p.loc[:, electric_ext]
            + p.loc[:, heat_ext]
            - p_nom.rename(rename).loc[electric_ext]
        )
        n.model.add_constraints(lhs <= 0, name="chplink-top_iso_fuel_line_ext")

    # top_iso_fuel_line for fixed
    if not electric_fix.empty:
        lhs = p.loc[:, electric_fix] + p.loc[:, heat_fix]
        rhs = n.links.p_nom[electric_fix]
        n.model.add_constraints(lhs <= rhs, name="chplink-top_iso_fuel_line_fix")

    # back-pressure
    if not electric.empty:
        lhs = (
            p.loc[:, heat] * (n.links.efficiency[heat] * n.links.c_b[electric].values)
            - p.loc[:, electric] * n.links.efficiency[electric]
        )
        n.model.add_constraints(lhs <= rhs, name="chplink-backpressure")


def add_pipe_retrofit_constraint(n):
    """
    Add constraint for retrofitting existing CH4 pipelines to H2 pipelines.
    """
    if "reversed" not in n.links.columns:
        n.links["reversed"] = False
    gas_pipes_i = n.links.query(
        "carrier == 'gas pipeline' and p_nom_extendable and ~reversed"
    ).index
    h2_retrofitted_i = n.links.query(
        "carrier == 'H2 pipeline retrofitted' and p_nom_extendable and ~reversed"
    ).index

    if h2_retrofitted_i.empty or gas_pipes_i.empty:
        return

    p_nom = n.model["Link-p_nom"]

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    lhs = p_nom.loc[gas_pipes_i] + CH4_per_H2 * p_nom.loc[h2_retrofitted_i]
    rhs = n.links.p_nom[gas_pipes_i]
    if not PYPSA_V1:
        rhs = rhs.rename_axis("Link-ext")

    n.model.add_constraints(lhs == rhs, name="Link-pipe_retrofit")


def add_flexible_egs_constraint(n):
    """
    Upper bounds the charging capacity of the geothermal reservoir according to
    the well capacity.
    """
    well_index = n.links.loc[n.links.carrier == "geothermal heat"].index
    storage_index = n.storage_units.loc[
        n.storage_units.carrier == "geothermal heat"
    ].index

    p_nom_rhs = n.model["Link-p_nom"].loc[well_index]
    p_nom_lhs = n.model["StorageUnit-p_nom"].loc[storage_index]

    n.model.add_constraints(
        p_nom_lhs <= p_nom_rhs,
        name="upper_bound_charging_capacity_of_geothermal_reservoir",
    )


def add_import_limit_constraint(n: pypsa.Network, sns: pd.DatetimeIndex):
    """
    Add constraint for limiting green energy imports (synthetic and biomass).
    Does not include fossil fuel imports.
    """

    nyears = n.snapshot_weightings.generators.sum() / 8760

    import_links = n.links.loc[n.links.carrier.str.contains("import")].index
    import_gens = n.generators.loc[n.generators.carrier.str.contains("import")].index

    limit = n.config["sector"]["imports"]["limit"]
    limit_sense = n.config["sector"]["imports"]["limit_sense"]

    if (import_links.empty and import_gens.empty) or not np.isfinite(limit):
        return

    weightings = n.snapshot_weightings.loc[sns, "generators"]

    # everything needs to be in MWh_fuel
    eff = n.links.loc[import_links, "efficiency"]

    p_gens = n.model["Generator-p"].loc[sns, import_gens]
    p_links = n.model["Link-p"].loc[sns, import_links]

    lhs = (p_gens * weightings).sum() + (p_links * eff * weightings).sum()

    rhs = limit * 1e6 * nyears

    n.model.add_constraints(lhs, limit_sense, rhs, name="import_limit")


def _get_diagnostic_slack_config(config: dict) -> dict:
    solving = config.get("solving", {})
    options = solving.get("options", {})
    slack_config = options.get("diagnostic_slacks", {})
    if not slack_config:
        slack_config = solving.get("diagnostic_slacks", {})

    if isinstance(slack_config, bool):
        slack_config = {"enable": slack_config}

    return slack_config or {}


def _diagnostic_slacks_enabled(config: dict) -> bool:
    return bool(_get_diagnostic_slack_config(config).get("enable", False))


def _compile_diagnostic_slack_patterns(slack_config: dict) -> list[re.Pattern]:
    patterns = list(slack_config.get("constraint_patterns", []))
    groups = slack_config.get("groups")

    if groups is None:
        groups = list(DEFAULT_DIAGNOSTIC_SLACK_PATTERNS)
    elif isinstance(groups, str):
        groups = [groups]

    for group in groups:
        patterns.extend(DEFAULT_DIAGNOSTIC_SLACK_PATTERNS.get(group, []))

    include_names = slack_config.get("constraint_names", [])
    if isinstance(include_names, str):
        include_names = [include_names]
    patterns.extend(f"^{re.escape(name)}$" for name in include_names)

    return [re.compile(pattern, re.IGNORECASE) for pattern in patterns]


def _constraint_coords(constraint: linopy.constraints.Constraint) -> dict:
    return {
        dim: constraint.coords[dim]
        for dim in constraint.sign.dims
        if dim in constraint.coords
    }


def _add_diagnostic_slack_variable(
    model: linopy.Model,
    name: str,
    coords: dict,
    mask: xr.DataArray | None,
) -> linopy.variables.Variable:
    if not coords:
        return model.add_variables(lower=0, name=name)

    coord_dims = set(coords)
    if mask is not None and not set(mask.dims).issubset(coord_dims):
        mask = None

    dims = list(coords)
    lower = xr.DataArray(0.0, coords=coords, dims=dims)
    return model.add_variables(
        lower=lower,
        mask=mask,
        name=name,
    )


def add_diagnostic_slacks(n: pypsa.Network) -> None:
    """
    Relax selected constraints with expensive non-negative slacks.

    This is a diagnostic tool for infeasible models. It should be enabled only
    for targeted runs via ``solving.options.diagnostic_slacks.enable``. The
    slacks are penalised in the objective and later exported so the binding
    source of infeasibility can be identified.
    """
    slack_config = _get_diagnostic_slack_config(n.config)
    if not slack_config.get("enable", False):
        return

    patterns = _compile_diagnostic_slack_patterns(slack_config)
    if not patterns:
        logger.warning("Diagnostic slacks enabled but no constraint patterns configured.")
        return

    penalty = float(slack_config.get("penalty", 1e9))
    model = n.model
    slack_specs = []
    objective_delta = None

    def add_to_objective(expr):
        nonlocal objective_delta
        objective_delta = expr if objective_delta is None else objective_delta + expr

    for name in list(model.constraints.data):
        if name.startswith("diagnostic_slack"):
            continue
        if not any(pattern.search(name) for pattern in patterns):
            continue

        constraint = model.constraints[name]
        sign_values = np.asarray(constraint.sign.values).reshape(-1)
        signs = pd.unique(pd.Series(sign_values).dropna().astype(str))
        if len(signs) != 1:
            logger.warning(
                "Skipping diagnostic slack for mixed-sign constraint block %s: %s",
                name,
                signs,
            )
            continue

        sign = signs[0]
        lhs = constraint.lhs
        rhs = constraint.rhs
        coords = _constraint_coords(constraint)
        mask = constraint.mask
        key = _sanitize_constraint_key(name)

        if sign in ("<=", "<", "≤"):
            slack_name = f"diagnostic_slack__{key}"
            slack = _add_diagnostic_slack_variable(
                model, name=slack_name, coords=coords, mask=mask
            )
            relaxed = lhs - slack <= rhs
            model.constraints.remove(name)
            model.add_constraints(relaxed, name=name)
            add_to_objective(penalty * slack.sum())
            slack_specs.append(
                {"constraint": name, "variable": slack_name, "direction": "upper"}
            )

        elif sign in (">=", ">", "≥"):
            slack_name = f"diagnostic_slack__{key}"
            slack = _add_diagnostic_slack_variable(
                model, name=slack_name, coords=coords, mask=mask
            )
            relaxed = lhs + slack >= rhs
            model.constraints.remove(name)
            model.add_constraints(relaxed, name=name)
            add_to_objective(penalty * slack.sum())
            slack_specs.append(
                {"constraint": name, "variable": slack_name, "direction": "lower"}
            )

        elif sign in ("=", "=="):
            slack_pos_name = f"diagnostic_slack_pos__{key}"
            slack_neg_name = f"diagnostic_slack_neg__{key}"
            slack_pos = _add_diagnostic_slack_variable(
                model, name=slack_pos_name, coords=coords, mask=mask
            )
            slack_neg = _add_diagnostic_slack_variable(
                model, name=slack_neg_name, coords=coords, mask=mask
            )
            relaxed = lhs + slack_pos - slack_neg == rhs
            model.constraints.remove(name)
            model.add_constraints(relaxed, name=name)
            add_to_objective(penalty * (slack_pos.sum() + slack_neg.sum()))
            slack_specs.extend(
                [
                    {
                        "constraint": name,
                        "variable": slack_pos_name,
                        "direction": "positive",
                    },
                    {
                        "constraint": name,
                        "variable": slack_neg_name,
                        "direction": "negative",
                    },
                ]
            )

        else:
            logger.warning(
                "Skipping diagnostic slack for constraint block %s with sign %s",
                name,
                sign,
            )

    if not slack_specs:
        logger.warning("Diagnostic slacks enabled but no matching constraints found.")
        return

    model.add_objective(
        model.objective.expression + objective_delta,
        overwrite=True,
        sense=model.objective.sense,
    )
    n._diagnostic_slack_specs = slack_specs
    logger.warning(
        "Added %s diagnostic slack variable blocks with objective penalty %.3g.",
        len(slack_specs),
        penalty,
    )


def export_diagnostic_slacks(
    n: pypsa.Network, solver_log: str | os.PathLike | None, tolerance: float = 1e-6
) -> None:
    specs = getattr(n, "_diagnostic_slack_specs", [])
    if not specs or solver_log is None:
        return

    solver_log = Path(solver_log)
    stem = solver_log.stem
    if stem.endswith("_solver"):
        stem = stem[: -len("_solver")]
    detail_path = solver_log.with_name(f"{stem}_diagnostic_slacks.csv")
    summary_path = solver_log.with_name(f"{stem}_diagnostic_slacks_summary.csv")

    rows = []
    summary_rows = []
    for spec in specs:
        try:
            solution = n.model[spec["variable"]].solution
            if solution.ndim == 0:
                frame = pd.DataFrame({"value": [float(solution.values)]})
            else:
                frame = solution.to_dataframe(name="value").reset_index()
        except Exception as exc:  # pragma: no cover - diagnostic fallback
            logger.warning(
                "Could not export diagnostic slack variable %s: %s",
                spec["variable"],
                exc,
            )
            continue

        if "value" not in frame:
            continue
        frame = frame.dropna(subset=["value"])
        nonzero = frame[frame["value"].abs() > tolerance].copy()

        summary_rows.append(
            {
                "constraint": spec["constraint"],
                "variable": spec["variable"],
                "direction": spec["direction"],
                "nonzero_count": int(len(nonzero)),
                "total_slack": float(nonzero["value"].abs().sum()),
                "max_slack": float(nonzero["value"].abs().max())
                if not nonzero.empty
                else 0.0,
            }
        )

        for _, row in nonzero.iterrows():
            coord_values = {
                column: row[column]
                for column in nonzero.columns
                if column != "value"
            }
            rows.append(
                {
                    "constraint": spec["constraint"],
                    "variable": spec["variable"],
                    "direction": spec["direction"],
                    "value": float(row["value"]),
                    **coord_values,
                }
            )

    pd.DataFrame(rows).to_csv(detail_path, index=False)
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary = summary.sort_values(["total_slack", "max_slack"], ascending=False)
    summary.to_csv(summary_path, index=False)
    logger.warning(
        "Wrote diagnostic slack reports to %s and %s", detail_path, summary_path
    )


def add_co2_atmosphere_constraint(n, snapshots):
    glcs = n.global_constraints[n.global_constraints.type == "co2_atmosphere"]

    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last_i = snapshots[-1]
            lhs = n.model["Store-e"].loc[last_i, stores.index]
            rhs = glc.constant

            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def extra_functionality(
    n: pypsa.Network, snapshots: pd.DatetimeIndex, planning_horizons: str | None = None
) -> None:
    """
    Add custom constraints and functionality.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance with config and params attributes
    snapshots : pd.DatetimeIndex
        Simulation timesteps
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments ``opts`` and
    ``snakemake.config`` are expected to be attached to the network.
    """
    config = n.config
    constraints = config["solving"].get("constraints", {})
    if constraints["BAU"] and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if constraints["SAFE"] and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if constraints["CCL"] and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config, planning_horizons)

    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)

    if EQ_o := constraints["EQ"]:
        add_EQ_constraints(n, EQ_o.replace("EQ", ""))

    if {"solar-hsat", "solar"}.issubset(
        config["electricity"]["renewable_carriers"]
    ) and {"solar-hsat", "solar"}.issubset(
        config["electricity"]["extendable_carriers"]["Generator"]
    ):
        add_solar_potential_constraints(n, config)

    limit_max_growth = config.get("sector", {}).get("limit_max_growth")
    if (
        not n._multi_invest
        and planning_horizons
        and limit_max_growth is not None
        and limit_max_growth.get("enable", False)
    ):
        add_max_growth_myopic(n, limit_max_growth, planning_horizons)

    if n.config.get("sector", {}).get("tes", False):
        if n.buses.index.str.contains(
            r"urban central heat|urban decentral heat|rural heat",
            case=False,
            na=False,
        ).any():
            add_TES_energy_to_power_ratio_constraints(n)
            add_TES_charger_ratio_constraints(n)

    add_battery_constraints(n)
    add_lossy_bidirectional_link_constraints(n)
    add_pipe_retrofit_constraint(n)
    if n._multi_invest:
        add_carbon_constraint(n, snapshots)
        add_carbon_budget_constraint(n, snapshots)
        add_retrofit_gas_boiler_constraint(n, snapshots)
    else:
        add_co2_atmosphere_constraint(n, snapshots)

    if config["sector"]["enhanced_geothermal"]["enable"]:
        add_flexible_egs_constraint(n)

    if config["sector"]["imports"]["enable"]:
        add_import_limit_constraint(n, snapshots)

    add_cdr_credit_accounting(n, planning_horizons)

    cdr_credit_limit = n.params.get("cdr_credit_limit")
    cdr_credit_limit_by_year = n.params.get("cdr_credit_limit_by_year")
    cdr_credit_scope = n.params.get("cdr_credit_scope") or []
    cdr_credit_timing = n.params.get("cdr_credit_timing", "capture")
    credit_limit = cdr_credit_limit_by_year or cdr_credit_limit
    if credit_limit and planning_horizons:
        add_cdr_credit_limit(
            n,
            limit_dict=credit_limit,
            planning_horizons=planning_horizons,
            cdr_credit_scope=cdr_credit_scope,
            cdr_credit_timing=cdr_credit_timing,
        )

    if n.params.custom_extra_functionality:
        source_path = n.params.custom_extra_functionality
        assert os.path.exists(source_path), f"{source_path} does not exist"
        sys.path.append(os.path.dirname(source_path))
        module_name = os.path.splitext(os.path.basename(source_path))[0]
        module = importlib.import_module(module_name)
        custom_extra_functionality = getattr(module, module_name)
        custom_extra_functionality(n, snapshots, snakemake)  # pylint: disable=E0601

    add_diagnostic_slacks(n)


def check_objective_value(n: pypsa.Network, solving: dict) -> None:
    """
    Check if objective value matches expected value within tolerance.

    Parameters
    ----------
    n : pypsa.Network
        Network with solved objective
    solving : Dict
        Dictionary containing objective checking parameters

    Raises
    ------
    ObjectiveValueError
        If objective value differs from expected value beyond tolerance
    """
    check_objective = solving["check_objective"]
    if check_objective["enable"]:
        atol = check_objective["atol"]
        rtol = check_objective["rtol"]
        expected_value = check_objective["expected_value"]
        if not np.isclose(n.objective, expected_value, atol=atol, rtol=rtol):
            raise ObjectiveValueError(
                f"Objective value {n.objective} differs from expected value "
                f"{expected_value} by more than {atol}."
            )


def collect_kwargs(
    config: dict,
    solving: dict,
    planning_horizons: str | None = None,
    log_fn: str | None = None,
    mode: str = "single",
) -> tuple[dict, dict]:
    """
    Prepare keyword arguments separated for model creation and model solving.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing solver settings
    solving : dict
        Dictionary of solving options and configuration
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    log_fn : str, optional
        Path to solver log file
    mode : str, optional
        Optimization mode: 'single', 'rolling_horizon', or 'iterative'
        Default is 'single'

    Returns
    -------
    tuple[dict, dict]
        Two dictionaries: (model_kwargs, solve_kwargs)
        - model_kwargs: Arguments for n.optimize.create_model()
        - solve_kwargs: Arguments for n.optimize.solve_model()
        For 'rolling_horizon' and 'iterative' modes, returns merged kwargs
        with additional mode-specific parameters
    """
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    # Model creation kwargs
    model_kwargs = {}
    model_kwargs["multi_investment_periods"] = config["foresight"] == "perfect"
    model_kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    model_kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment", False
    )

    # Solve kwargs
    solver_name = solving["solver"]["name"]
    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}

    io_api = cf_solving.get("io_api", None)
    if io_api is None and solver_name == "gurobi":
        logger.info("No io_api configured; defaulting to Gurobi MPS interface.")
        io_api = "mps"

    if io_api == "mps":
        patch_linopy_highspy_name_export()

    solve_kwargs = {}
    solve_kwargs["solver_name"] = solver_name
    solve_kwargs["solver_options"] = solver_options
    solve_kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)
    solve_kwargs["io_api"] = io_api
    solve_kwargs["keep_files"] = cf_solving.get("keep_files", False)

    if log_fn:
        solve_kwargs["log_fn"] = log_fn

    oetc = solving.get("oetc", None)
    if oetc:
        oetc["credentials"] = OetcCredentials(
            email=os.environ["OETC_EMAIL"], password=os.environ["OETC_PASSWORD"]
        )
        oetc["solver"] = solver_name
        oetc["solver_options"] = solver_options
        oetc_settings = OetcSettings(**oetc)
        oetc_handler = OetcHandler(oetc_settings)
        solve_kwargs["remote"] = oetc_handler

    if solver_name == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    # Handle special modes
    if mode == "rolling_horizon":
        all_kwargs = {**model_kwargs, **solve_kwargs}
        all_kwargs["horizon"] = cf_solving.get("horizon", 365)
        all_kwargs["overlap"] = cf_solving.get("overlap", 0)
        return all_kwargs, {}

    elif mode == "iterative":
        all_kwargs = {**model_kwargs, **solve_kwargs}
        all_kwargs["track_iterations"] = cf_solving["track_iterations"]
        all_kwargs["min_iterations"] = cf_solving["min_iterations"]
        all_kwargs["max_iterations"] = cf_solving["max_iterations"]

        if cf_solving["post_discretization"].get("enable", False):
            logger.info("Add post-discretization parameters.")
            all_kwargs.update(cf_solving["post_discretization"])

        return all_kwargs, {}

    return model_kwargs, solve_kwargs


def create_optimization_model(
    n: pypsa.Network,
    config: dict,
    params: dict,
    model_kwargs: dict,
    solve_kwargs: dict,
    planning_horizons: str | None = None,
) -> None:
    """
    Prepare optimization problem by creating model and adding extra functionality.

    This function:
    1. Attaches config and params to network for extra_functionality
    2. Creates the optimization model
    3. Adds extra functionality (custom constraints)

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : dict
        Configuration dictionary containing solver settings
    params : dict
        Dictionary of solving parameters
    model_kwargs : dict
        Arguments for n.optimize.create_model()
    solve_kwargs : dict
        Arguments for n.optimize.solve_model()
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    """
    # Add config and params to network for extra_functionality
    n.config = config
    n.params = params

    # Create optimization model
    logger.info("Creating optimization model...")
    n.optimize.create_model(**model_kwargs)

    # Add extra functionality (custom constraints)
    logger.info("Adding extra functionality (custom constraints)...")
    extra_functionality(n, n.snapshots, planning_horizons)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_sector_network",
            opts="",
            clusters="5",
            configfiles="config/test/config.overnight.yaml",
            sector_opts="",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]
    cf_solving = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    # Load network
    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    # Prepare network (settings before solving)
    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
    )

    # Determine solve mode
    rolling_horizon = cf_solving.get("rolling_horizon", False)
    skip_iterations = cf_solving.get("skip_iterations", False)

    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )

    # Solve network based on mode
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        if rolling_horizon and snakemake.rule == "solve_operations_network":
            logger.info("Using rolling horizon optimization...")
            all_kwargs, _ = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="rolling_horizon",
            )

            n.config = snakemake.config
            n.params = snakemake.params
            all_kwargs["extra_functionality"] = partial(
                extra_functionality, planning_horizons=planning_horizons
            )
            n.optimize.optimize_with_rolling_horizon(**all_kwargs)
            status, condition = "", ""

        elif skip_iterations:
            logger.info("Using single-pass optimization...")
            model_kwargs, solve_kwargs = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="single",
            )
            create_optimization_model(
                n,
                config=snakemake.config,
                params=snakemake.params,
                model_kwargs=model_kwargs,
                solve_kwargs=solve_kwargs,
                planning_horizons=planning_horizons,
            )

            logger.info("Solving model...")
            status, condition = n.optimize.solve_model(**solve_kwargs)

        else:
            logger.info("Using iterative transmission expansion optimization...")

            all_kwargs, _ = collect_kwargs(
                snakemake.config,
                snakemake.params.solving,
                planning_horizons,
                log_fn=snakemake.log.solver,
                mode="iterative",
            )

            n.config = snakemake.config
            n.params = snakemake.params
            all_kwargs["extra_functionality"] = partial(
                extra_functionality, planning_horizons=planning_horizons
            )
            status, condition = n.optimize.optimize_transmission_expansion_iteratively(
                **all_kwargs
            )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    if _diagnostic_slacks_enabled(snakemake.config):
        export_diagnostic_slacks(
            n,
            solver_log=getattr(snakemake.log, "solver", None),
            tolerance=float(
                _get_diagnostic_slack_config(snakemake.config).get("tolerance", 1e-6)
            ),
        )
        raise RuntimeError(
            "Diagnostic slack run completed; not exporting relaxed solution as "
            "a solved network."
        )

    # Check results
    condition_str = str(condition).lower()
    if not rolling_horizon:
        if status != "ok":
            logger.warning(
                f"Solving status '{status}' with termination condition '{condition}'"
            )
        check_objective_value(n, snakemake.params.solving)

    if "infeasible" in condition_str:
        labels = n.model.compute_infeasibilities()
        logger.info(f"Labels:\n{labels}")
        n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'. Infeasibilities computed.")

    if status != "ok":
        raise RuntimeError(
            f"Solving failed with status '{status}' and termination condition "
            f"'{condition}'. Discarding solution and skipping exports."
        )

    if "warning" in condition_str:
        raise RuntimeError("Solving status 'warning'. Discarding solution.")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    with open(snakemake.output.config, "w") as file:
        yaml.dump(
            n.meta,
            file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )

    if hasattr(snakemake.output, "cdr_credit_accounting"):
        try:
            export_cdr_credit_accounting(
                n,
                output_path=snakemake.output.cdr_credit_accounting,
                planning_horizons=planning_horizons,
            )
        except Exception as exc:
            logger.exception("Failed to export CDR credit accounting.")
            export_cdr_credit_accounting_failure(
                n,
                output_path=snakemake.output.cdr_credit_accounting,
                planning_horizons=planning_horizons,
                error=exc,
            )
