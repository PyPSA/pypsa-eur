# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build a stochastic PyPSA network from a deterministic pre-solve network.

This script is intended to run before solve_network.py:
- load the deterministic pre-solve network;
- apply prepare_network(...);
- call n.set_scenarios(...);
- apply scenario-specific modifications;
- export the stochastic pre-solve network.
"""

import logging
import re
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import pypsa
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts._helpers import (  # noqa: E402
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import prepare_network  # noqa: E402


logger = logging.getLogger(__name__)


# ---------------------------
# YAML / generic helpers
# ---------------------------

def _read_yaml_maybe(path: str | None) -> dict:
    """Read a YAML file if path is provided and exists; return {} otherwise."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Stochastic config file not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _ensure_dict(x: Any, name: str) -> dict:
    """Ensure x is a dict, else raise."""
    if x is None:
        return {}
    if isinstance(x, dict):
        return x
    raise TypeError(f"{name} must be a dict; got {type(x).__name__}")


def _get_level_names(idx: pd.Index) -> pd.Index:
    """Return the 'name' level if MultiIndex, else the index itself."""
    if isinstance(idx, pd.MultiIndex):
        if "name" in idx.names:
            return idx.get_level_values("name")
        return idx.get_level_values(-1)
    return idx


def _merge_stochastic_param(stochastic_param: dict) -> dict:
    """Merge inline stochastic config with optional external YAML."""
    stoch = _ensure_dict(stochastic_param, "stochastic_scenarios")
    external = _read_yaml_maybe(stoch.get("file", ""))
    if external:
        merged = dict(external)
        merged.update(stoch)
        stoch = merged
    return stoch


# ---------------------------
# Low-level helpers for loads/links
# ---------------------------

def _base_component_table(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of a component table with scenario removed from the index if present."""
    base = df.copy()
    base.index = _get_level_names(base.index)
    return base


def _base_loads_table(n: pypsa.Network) -> pd.DataFrame:
    """Return loads table indexed only by load name."""
    return _base_component_table(n.loads)


def _base_links_table(n: pypsa.Network) -> pd.DataFrame:
    """Return links table indexed only by link name."""
    return _base_component_table(n.links)


def _ts_column_key(ts: pd.DataFrame, name: str, scenario: str | None = None):
    """
    Return the correct column key for a time series table.

    For deterministic tables the key is the plain component name.
    For stochastic tables the key is the tuple (scenario, name).
    """
    if not isinstance(ts.columns, pd.MultiIndex):
        if name not in ts.columns:
            raise KeyError(f"Column '{name}' not found in deterministic time series table.")
        return name

    if scenario is None:
        raise ValueError(
            f"Time series table is stochastic but no scenario was provided for column '{name}'."
        )

    if "scenario" in ts.columns.names:
        scenarios = ts.columns.get_level_values("scenario")
    else:
        scenarios = ts.columns.get_level_values(0)

    if "name" in ts.columns.names:
        names = ts.columns.get_level_values("name")
    else:
        names = ts.columns.get_level_values(1)

    mask = (scenarios == scenario) & (names == name)
    matches = ts.columns[mask]
    if len(matches) == 0:
        raise KeyError(f"Column for scenario='{scenario}', name='{name}' not found.")
    if len(matches) > 1:
        raise ValueError(f"Multiple columns found for scenario='{scenario}', name='{name}'.")
    return matches[0]


def _read_ts_series(ts: pd.DataFrame, name: str, scenario: str | None = None) -> pd.Series:
    """Return a copy of one time series column."""
    key = _ts_column_key(ts, name=name, scenario=scenario)
    return ts.loc[:, key].copy()


def _write_ts_series(
    ts: pd.DataFrame,
    name: str,
    values: pd.Series | np.ndarray,
    scenario: str | None = None,
) -> None:
    """Overwrite one time series column."""
    key = _ts_column_key(ts, name=name, scenario=scenario)

    if isinstance(values, pd.Series):
        v = values.reindex(ts.index)
        if v.isnull().any():
            raise ValueError(f"NaNs encountered while writing series '{name}'.")
        ts.loc[:, key] = v.values
    else:
        arr = np.asarray(values)
        if arr.ndim != 1 or len(arr) != len(ts.index):
            raise ValueError(
                f"Invalid shape for series '{name}': expected ({len(ts.index)},), got {arr.shape}."
            )
        ts.loc[:, key] = arr

def _scale_loads_by_carrier(
    n: pypsa.Network,
    carrier: str,
    factor: float,
    scenario: str | None = None,
) -> None:
    """Scale all loads belonging to a given carrier."""
    names = _find_load_names_by_carrier(n, carrier)
    if not names:
        raise KeyError(f"No loads found for carrier '{carrier}'.")

    for name in names:
        s = _read_load_series(n, name, scenario=scenario)
        _write_load_series(n, name, s * factor, scenario=scenario)

    logger.info(
        "Scaled %s load(s) with carrier '%s' by factor %.4f.",
        len(names),
        carrier,
        factor,
    )

def _read_link_efficiency_series(
    n: pypsa.Network,
    link_name: str,
    scenario: str | None = None,
) -> pd.Series:
    """Return one link efficiency time series."""
    return _read_ts_series(n.links_t.efficiency, name=link_name, scenario=scenario)


def _carrier_names_from_table(df: pd.DataFrame, carrier: str) -> list[str]:
    """Return component names matching a carrier from a base component table."""
    base = _base_component_table(df)
    if "carrier" not in base.columns:
        return []
    mask = base["carrier"].astype(str).eq(carrier)
    return base.index[mask].tolist()


def _find_load_names_by_carrier(n: pypsa.Network, carrier: str) -> list[str]:
    """Return load names whose carrier matches the requested value."""
    return _carrier_names_from_table(n.loads, carrier)


def _find_link_names_by_carrier(n: pypsa.Network, carrier: str) -> list[str]:
    """Return link names whose carrier matches the requested value."""
    return _carrier_names_from_table(n.links, carrier)


def _extract_prefix(name: str, suffix: str) -> str:
    """Strip a known suffix from a component name and return the prefix."""
    if not name.endswith(suffix):
        raise ValueError(f"Name '{name}' does not end with suffix '{suffix}'.")
    return name[: -len(suffix)]


def _assert_load_exists(n: pypsa.Network, load_name: str) -> None:
    """Raise if a load is missing."""
    loads = _base_loads_table(n)
    if load_name not in loads.index:
        raise KeyError(f"Required load '{load_name}' not found in n.loads.")


def _assert_link_exists(n: pypsa.Network, link_name: str) -> None:
    """Raise if a link is missing."""
    links = _base_links_table(n)
    if link_name not in links.index:
        raise KeyError(f"Required link '{link_name}' not found in n.links.")

def _load_row_key(n: pypsa.Network, load_name: str, scenario: str | None = None):
    """
    Return the correct row key for n.loads.

    For deterministic tables the key is the plain load name.
    For stochastic static tables the key is the tuple (scenario, name).
    """
    if not isinstance(n.loads.index, pd.MultiIndex):
        if load_name not in n.loads.index:
            raise KeyError(f"Load '{load_name}' not found in deterministic n.loads.")
        return load_name

    if scenario is None:
        raise ValueError(
            f"n.loads has a stochastic MultiIndex but no scenario was provided for '{load_name}'."
        )

    if "scenario" in n.loads.index.names:
        scenarios = n.loads.index.get_level_values("scenario")
    else:
        scenarios = n.loads.index.get_level_values(0)

    if "name" in n.loads.index.names:
        names = n.loads.index.get_level_values("name")
    else:
        names = n.loads.index.get_level_values(1)

    mask = (scenarios == scenario) & (names == load_name)
    matches = n.loads.index[mask]
    if len(matches) == 0:
        raise KeyError(f"Static load row for scenario='{scenario}', name='{load_name}' not found.")
    if len(matches) > 1:
        raise ValueError(
            f"Multiple static load rows found for scenario='{scenario}', name='{load_name}'."
        )
    return matches[0]


def _has_load_timeseries_column(
    n: pypsa.Network,
    load_name: str,
    scenario: str | None = None,
) -> bool:
    """Return True if the load exists in n.loads_t.p_set."""
    ts = n.loads_t.p_set
    if not isinstance(ts.columns, pd.MultiIndex):
        return load_name in ts.columns

    if scenario is None:
        return False

    if "scenario" in ts.columns.names:
        scenarios = ts.columns.get_level_values("scenario")
    else:
        scenarios = ts.columns.get_level_values(0)

    if "name" in ts.columns.names:
        names = ts.columns.get_level_values("name")
    else:
        names = ts.columns.get_level_values(1)

    return ((scenarios == scenario) & (names == load_name)).any()


def _read_load_series(n: pypsa.Network, load_name: str, scenario: str | None = None) -> pd.Series:
    """
    Return the effective load time series.

    Priority:
    1. n.loads_t.p_set column if present
    2. broadcast static n.loads.p_set over all snapshots
    """
    if _has_load_timeseries_column(n, load_name, scenario=scenario):
        return _read_ts_series(n.loads_t.p_set, name=load_name, scenario=scenario)

    row_key = _load_row_key(n, load_name, scenario=scenario)
    value = n.loads.loc[row_key, "p_set"]
    return pd.Series(float(value), index=n.snapshots)


def _write_static_load_value(
    n: pypsa.Network,
    load_name: str,
    value: float,
    scenario: str | None = None,
) -> None:
    """Write a scalar value into n.loads.p_set."""
    row_key = _load_row_key(n, load_name, scenario=scenario)
    n.loads.loc[row_key, "p_set"] = float(value)


def _ensure_load_timeseries_column(
    n: pypsa.Network,
    load_name: str,
    scenario: str | None = None,
) -> None:
    """
    Ensure that n.loads_t.p_set contains a column for the requested load.

    If missing, initialize it from the current effective series.
    """
    if _has_load_timeseries_column(n, load_name, scenario=scenario):
        return

    base_series = _read_load_series(n, load_name, scenario=scenario)

    ts = n.loads_t.p_set
    if not isinstance(ts.columns, pd.MultiIndex):
        ts.loc[:, load_name] = base_series.values
        ts.columns.name = "name"
        return

    if scenario is None:
        raise ValueError(
            f"Cannot create stochastic time series column for '{load_name}' without scenario."
        )

    new_col = pd.MultiIndex.from_tuples(
        [(scenario, load_name)],
        names=ts.columns.names,
    )
    new_df = pd.DataFrame(base_series.values, index=ts.index, columns=new_col)
    n.loads_t.p_set = pd.concat([ts, new_df], axis=1)


def _write_load_series(
    n: pypsa.Network,
    load_name: str,
    values: pd.Series | np.ndarray,
    scenario: str | None = None,
) -> None:
    """
    Write a load profile.

    - If the profile is constant and the load has no existing timeseries column, write to static n.loads.p_set
    - Otherwise create/use a column in n.loads_t.p_set
    """
    if not isinstance(values, pd.Series):
        arr = np.asarray(values)
        if arr.ndim != 1 or len(arr) != len(n.snapshots):
            raise ValueError(
                f"Invalid shape for series '{load_name}': expected ({len(n.snapshots)},), got {arr.shape}."
            )
        values = pd.Series(arr, index=n.snapshots)
    else:
        values = values.reindex(n.snapshots)
        if values.isnull().any():
            raise ValueError(f"NaNs encountered while writing series '{load_name}'.")

    is_constant = np.allclose(values.values, values.values[0])

    if is_constant and not _has_load_timeseries_column(n, load_name, scenario=scenario):
        _write_static_load_value(n, load_name, float(values.iloc[0]), scenario=scenario)
        return

    _ensure_load_timeseries_column(n, load_name, scenario=scenario)
    _write_ts_series(n.loads_t.p_set, name=load_name, values=values, scenario=scenario)

# ---------------------------
# Structured scenario builders
# ---------------------------

TRANSPORT_ELECTRIC_EFFICIENCY = 53.19
TRANSPORT_ICE_EFFICIENCY = 16.0712
URBAN_HEAT_CENTRAL_ALPHA = 0.98


def _scenario_agriculture_full_electric(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """
    Electrify agriculture:
    - move all agriculture machinery oil demand to agriculture machinery electric (1:1)
    - convert agriculture heat demand into agriculture electricity using local rural air heat pump COP
    """
    del config  # unused for now

    # Part 1: machinery oil -> machinery electric
    oil_names = _find_load_names_by_carrier(n, "agriculture machinery oil")
    if not oil_names:
        raise KeyError("No loads found for carrier 'agriculture machinery oil'.")

    count_machinery = 0
    for oil_name in oil_names:
        prefix = _extract_prefix(oil_name, "agriculture machinery oil")
        elec_name = f"{prefix}agriculture machinery electric"

        _assert_load_exists(n, elec_name)

        oil_series = _read_load_series(n, oil_name, scenario=scenario)
        elec_series = _read_load_series(n, elec_name, scenario=scenario)

        _write_load_series(n, elec_name, elec_series + oil_series, scenario=scenario)
        _write_load_series(
            n, oil_name, pd.Series(0.0, index=n.loads_t.p_set.index), scenario=scenario
        )
        count_machinery += 1

    # Part 2: agriculture heat -> agriculture electricity via rural air heat pump COP
    heat_names = _find_load_names_by_carrier(n, "agriculture heat")
    if not heat_names:
        raise KeyError("No loads found for carrier 'agriculture heat'.")

    count_heat = 0
    for heat_name in heat_names:
        prefix = _extract_prefix(heat_name, "agriculture heat")
        elec_name = f"{prefix}agriculture electricity"
        hp_name = f"{prefix}rural air heat pump"

        _assert_load_exists(n, elec_name)
        _assert_link_exists(n, hp_name)

        heat_series = _read_load_series(n, heat_name, scenario=scenario)
        elec_series = _read_load_series(n, elec_name, scenario=scenario)
        cop = _read_link_efficiency_series(n, hp_name, scenario=scenario)

        if (cop <= 0).any():
            raise ValueError(f"Non-positive COP detected for link '{hp_name}'.")

        added_electricity = heat_series / cop

        _write_load_series(n, elec_name, elec_series + added_electricity, scenario=scenario)
        _write_load_series(
            n, heat_name, pd.Series(0.0, index=n.loads_t.p_set.index), scenario=scenario
        )
        count_heat += 1

    logger.info(
        "Applied scenario 'agriculture_full_electric' to %s machinery node(s) and %s heat node(s)%s.",
        count_machinery,
        count_heat,
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_agriculture_machinery_full_oil(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """Move all agriculture machinery electric demand to agriculture machinery oil (1:1)."""
    del config  # unused for now

    elec_names = _find_load_names_by_carrier(n, "agriculture machinery electric")
    if not elec_names:
        raise KeyError("No loads found for carrier 'agriculture machinery electric'.")

    count = 0
    for elec_name in elec_names:
        prefix = _extract_prefix(elec_name, "agriculture machinery electric")
        oil_name = f"{prefix}agriculture machinery oil"

        _assert_load_exists(n, oil_name)

        elec_series = _read_load_series(n, elec_name, scenario=scenario)
        oil_series = _read_load_series(n, oil_name, scenario=scenario)

        _write_load_series(n, oil_name, oil_series + elec_series, scenario=scenario)
        _write_load_series(
            n, elec_name, pd.Series(0.0, index=n.loads_t.p_set.index), scenario=scenario
        )
        count += 1

    logger.info(
        "Applied scenario 'agriculture_machinery_full_oil' to %s node(s)%s.",
        count,
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_shipping_full_methanol(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """
    Move all shipping oil demand to the global EU shipping methanol load (1:1).

    Current PyPSA-Eur structure:
    - shipping oil is nodal
    - shipping methanol is represented by a single global load: 'EU shipping methanol'
    """
    del config  # unused for now

    oil_names = _find_load_names_by_carrier(n, "shipping oil")
    if not oil_names:
        raise KeyError("No loads found for carrier 'shipping oil'.")

    methanol_name = "EU shipping methanol"
    _assert_load_exists(n, methanol_name)

    methanol_series = _read_load_series(n, methanol_name, scenario=scenario)
    total_oil = pd.Series(0.0, index=n.snapshots)

    for oil_name in oil_names:
        total_oil = total_oil + _read_load_series(n, oil_name, scenario=scenario)

    _write_load_series(n, methanol_name, methanol_series + total_oil, scenario=scenario)

    zero = pd.Series(0.0, index=n.snapshots)
    for oil_name in oil_names:
        _write_load_series(n, oil_name, zero, scenario=scenario)

    logger.info(
        "Applied scenario 'shipping_full_methanol': moved %s nodal shipping oil load(s) into '%s'%s.",
        len(oil_names),
        methanol_name,
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_urban_heat_full_central(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """Shift 98% of urban decentral heat demand to urban central heat."""
    del config  # unused for now

    decentral_names = _find_load_names_by_carrier(n, "urban decentral heat")
    if not decentral_names:
        raise KeyError("No loads found for carrier 'urban decentral heat'.")

    count = 0
    for decentral_name in decentral_names:
        prefix = _extract_prefix(decentral_name, "urban decentral heat")
        central_name = f"{prefix}urban central heat"

        _assert_load_exists(n, central_name)

        decentral_series = _read_load_series(n, decentral_name, scenario=scenario)
        central_series = _read_load_series(n, central_name, scenario=scenario)

        moved = URBAN_HEAT_CENTRAL_ALPHA * decentral_series
        remaining = (1.0 - URBAN_HEAT_CENTRAL_ALPHA) * decentral_series

        _write_load_series(n, central_name, central_series + moved, scenario=scenario)
        _write_load_series(n, decentral_name, remaining, scenario=scenario)
        count += 1

    logger.info(
        "Applied scenario 'urban_heat_full_central' with alpha=%.2f to %s node(s)%s.",
        URBAN_HEAT_CENTRAL_ALPHA,
        count,
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_land_transport_linear_ev(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """
    Reallocate land transport useful service:
    - 60% EV
    - 40% oil / ICE
    preserving useful transport service locally at each node.
    """
    del config  # unused for now

    ev_names = _find_load_names_by_carrier(n, "land transport EV")
    if not ev_names:
        raise KeyError("No loads found for carrier 'land transport EV'.")

    count = 0
    for ev_name in ev_names:
        prefix = _extract_prefix(ev_name, "land transport EV")
        oil_name = f"{prefix}land transport oil"

        _assert_load_exists(n, oil_name)

        ev_series = _read_load_series(n, ev_name, scenario=scenario)
        oil_series = _read_load_series(n, oil_name, scenario=scenario)

        useful_service = (
            ev_series * TRANSPORT_ELECTRIC_EFFICIENCY
            + oil_series * TRANSPORT_ICE_EFFICIENCY
        )

        ev_new = 0.60 * useful_service / TRANSPORT_ELECTRIC_EFFICIENCY
        oil_new = 0.40 * useful_service / TRANSPORT_ICE_EFFICIENCY

        _write_load_series(n, ev_name, ev_new, scenario=scenario)
        _write_load_series(n, oil_name, oil_new, scenario=scenario)
        count += 1

    logger.info(
        "Applied scenario 'land_transport_linear_ev' to %s node(s)%s.",
        count,
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_electricity_optimistic(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """Increase generic electricity demand by 10%."""
    del config  # unused for now

    _scale_loads_by_carrier(n, carrier="electricity", factor=1.10, scenario=scenario)

    logger.info(
        "Applied scenario 'electricity_optimistic'%s.",
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_industry_h2(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """
    Move gas for industry and solid biomass for industry demand to H2 for industry (1:1).
    """
    del config  # unused for now

    gas_names = _find_load_names_by_carrier(n, "gas for industry")
    if not gas_names:
        raise KeyError("No loads found for carrier 'gas for industry'.")

    count = 0
    for gas_name in gas_names:
        prefix = _extract_prefix(gas_name, "gas for industry")
        biomass_name = f"{prefix}solid biomass for industry"
        h2_name = f"{prefix}H2 for industry"

        _assert_load_exists(n, biomass_name)
        _assert_load_exists(n, h2_name)

        gas_series = _read_load_series(n, gas_name, scenario=scenario)
        biomass_series = _read_load_series(n, biomass_name, scenario=scenario)
        h2_series = _read_load_series(n, h2_name, scenario=scenario)

        _write_load_series(
            n,
            h2_name,
            h2_series + gas_series + biomass_series,
            scenario=scenario,
        )
        _write_load_series(
            n, gas_name, pd.Series(0.0, index=n.loads_t.p_set.index), scenario=scenario
        )
        _write_load_series(
            n, biomass_name, pd.Series(0.0, index=n.loads_t.p_set.index), scenario=scenario
        )
        count += 1

    logger.info(
        "Applied scenario 'industry_h2' to %s node(s)%s.",
        count,
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )

def _scenario_base(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """Base scenario with no modifications."""
    del n, config  # unused for now

    logger.info(
        "Applied scenario 'base' with no modifications%s.",
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _apply_named_structured_scenario(
    n: pypsa.Network,
    scenario_name: str,
    config: dict | None,
    scenario: str | None = None,
) -> None:
    """Dispatch a structured scenario by name."""
    if scenario_name not in STRUCTURED_SCENARIOS:
        known = ", ".join(sorted(STRUCTURED_SCENARIOS))
        raise ValueError(
            f"Unknown structured scenario '{scenario_name}'. Known structured scenarios: {known}"
        )

    logger.info(
        "Applying structured scenario builder '%s' (target=%s).",
        scenario_name,
        f"stochastic:{scenario}" if scenario is not None else "deterministic",
    )
    STRUCTURED_SCENARIOS[scenario_name](n=n, scenario=scenario, config=config or {})


def _resolve_deterministic_structured_scenario(
    config: dict,
    wildcards: Mapping[str, Any] | None = None,
) -> tuple[str | None, dict]:
    """
    Resolve the active structured scenario in deterministic mode.

    Precedence:
    1. config['structured_scenario']
    2. config['scenario']['structured_name']
    3. wildcards['run']
    4. config['run']['name']

    Returns
    -------
    scenario_name : str | None
    scenario_params : dict
    """
    wildcards = wildcards or {}

    scenario_block = config.get("scenario", {})
    if not isinstance(scenario_block, dict):
        scenario_block = {}

    run_block = config.get("run", {})
    if not isinstance(run_block, dict):
        run_block = {}

    candidates = [
        config.get("structured_scenario"),
        scenario_block.get("structured_name"),
        wildcards.get("run"),
        run_block.get("name"),
    ]

    for candidate in candidates:
        try:
            scenario_name, scenario_params = _normalize_structured_scenario_spec(candidate)
        except (TypeError, ValueError):
            continue

        if isinstance(scenario_name, str) and scenario_name in STRUCTURED_SCENARIOS:
            return scenario_name, scenario_params

    return None, {}


def _validate_stochastic_structured_scenarios(scenarios: Mapping[str, Any]) -> None:
    """
    Ensure all stochastic structured_scenario specifications are valid.

    Supported stochastic scenario formats
    ------------------------------------
    Old format:
        scenarios:
          base: 0.125
          scenario_a: 0.125

    New format:
        scenarios:
          base:
            prob: 0.125
            structured_scenario: null
          gas_expensive:
            prob: 0.125
            structured_scenario:
              name: modify_components
              params: {...}
    """
    for sc_name, sc_spec in scenarios.items():
        if np.isscalar(sc_spec):
            # Legacy format: scenario name itself is used as structured scenario name
            structured_spec = sc_name
        elif isinstance(sc_spec, dict):
            structured_spec = sc_spec.get("structured_scenario", sc_name)
        else:
            raise TypeError(
                f"Invalid stochastic scenario specification for '{sc_name}': "
                f"expected scalar probability or dict, got {type(sc_spec).__name__}"
            )

        scenario_name, _ = _normalize_structured_scenario_spec(structured_spec)

        if scenario_name is None:
            continue

        if scenario_name not in STRUCTURED_SCENARIOS:
            known = ", ".join(sorted(STRUCTURED_SCENARIOS))
            raise ValueError(
                f"Unknown structured scenario '{scenario_name}' declared for stochastic "
                f"scenario '{sc_name}'. Known structured scenarios: {known}"
            )

def _normalize_stochastic_scenarios_definition(
    scenarios: Mapping[str, Any],
) -> tuple[dict[str, float], dict[str, tuple[str | None, dict]]]:
    """
    Normalize stochastic scenario definitions.

    Returns
    -------
    probabilities : dict[str, float]
        Scenario probabilities for n.set_scenarios(...)
    structured_specs : dict[str, tuple[str | None, dict]]
        Mapping scenario_name -> (structured_scenario_name, params)

    Notes
    -----
    Supported formats:

    1. Legacy
       scenarios:
         base: 0.125
         agriculture_full_electric: 0.125

    2. Extended
       scenarios:
         base:
           prob: 0.125
           structured_scenario: null
         custom_case:
           prob: 0.125
           structured_scenario:
             name: modify_components
             params:
               rules: [...]
    """
    probabilities = {}
    structured_specs = {}

    for sc_name, sc_spec in scenarios.items():
        if np.isscalar(sc_spec):
            probabilities[sc_name] = float(sc_spec)
            structured_specs[sc_name] = (sc_name, {})
            continue

        if not isinstance(sc_spec, dict):
            raise TypeError(
                f"Invalid stochastic scenario specification for '{sc_name}': "
                f"expected scalar or dict, got {type(sc_spec).__name__}"
            )

        if "prob" not in sc_spec:
            raise ValueError(
                f"Stochastic scenario '{sc_name}' is missing required key 'prob'."
            )

        probabilities[sc_name] = float(sc_spec["prob"])

        structured_spec = sc_spec.get("structured_scenario", sc_name)
        structured_specs[sc_name] = _normalize_structured_scenario_spec(structured_spec)

    return probabilities, structured_specs

def _apply_structured_scenarios(
    n: pypsa.Network,
    config: dict,
    stochastic_param: dict,
    wildcards: Mapping[str, Any] | None = None,
) -> tuple[bool, list[str]]:
    """
    Apply structured scenarios in stochastic or deterministic mode.

    Returns
    -------
    is_stochastic : bool
        Whether the stochastic mode was enabled.
    active_names : list[str]
        List of structured scenario names that were applied.
    """
    wildcards = wildcards or {}
    stoch = _merge_stochastic_param(stochastic_param)
    enabled = bool(stoch.get("enable", stoch.get("enabled", False)))

    if enabled:
        scenarios = stoch.get("scenarios", None)
        if scenarios is None:
            raise ValueError(
                "stochastic_scenarios.enable=true but no scenarios were provided "
                "(inline or through the referenced YAML file)."
            )

        _validate_stochastic_structured_scenarios(scenarios)
        probabilities, structured_specs = _normalize_stochastic_scenarios_definition(scenarios)

        logger.info("Enabling stochastic scenarios via n.set_scenarios(...)")
        n.set_scenarios(probabilities)

        active_names = []
        logger.info("Structured stochastic scenarios detected: %s", list(probabilities.keys()))

        for sc in probabilities:
            scenario_name, scenario_params = structured_specs[sc]

            if scenario_name is None:
                logger.info(
                    "No structured scenario declared for stochastic scenario '%s'; skipping builder.",
                    sc,
                )
                continue

            _apply_named_structured_scenario(
                n=n,
                scenario_name=scenario_name,
                config=scenario_params,
                scenario=sc,
            )
            active_names.append(f"{sc}->{scenario_name}")

        return True, active_names

    scenario_name, scenario_params = _resolve_deterministic_structured_scenario(
        config=config,
        wildcards=wildcards,
    )

    if scenario_name is None:
        logger.info(
            "Stochastic mode disabled and no deterministic structured scenario detected. "
            "No structured scenario builder will be applied."
        )
        return False, []

    logger.info("Deterministic structured scenario detected: %s", scenario_name)
    _apply_named_structured_scenario(
        n=n,
        scenario_name=scenario_name,
        config=scenario_params,
        scenario=None,
    )
    return False, [scenario_name]


# ---------------------------
# Patch-based selectors/helpers
# ---------------------------

def _select_names_from_component(
    n: pypsa.Network,
    comp: str,
    selector: Mapping[str, Any],
) -> list[str]:
    """
    Select component names (without scenario level) based on a selector.

    Supported selector keys:
    - names: list[str] or str (regex)
    - carrier: str or list[str]
    - bus / bus0 / bus1: str or list[str] (exact match)
    - any other column in the component table: exact match
    """
    selector = _ensure_dict(selector, "selector")
    df = getattr(n, comp)

    idx_names = _get_level_names(df.index)
    base = df.copy()
    base.index = idx_names

    mask = pd.Series(True, index=base.index)

    names_sel = selector.get("names")
    if isinstance(names_sel, str):
        pattern = re.compile(names_sel)
        mask &= base.index.to_series().apply(lambda s: bool(pattern.search(s)))
    elif isinstance(names_sel, (list, tuple, set)):
        mask &= base.index.isin(list(names_sel))

    carrier_sel = selector.get("carrier")
    if carrier_sel is not None and "carrier" in base.columns:
        if isinstance(carrier_sel, str):
            carrier_sel = [carrier_sel]
        mask &= base["carrier"].isin(list(carrier_sel))

    for bcol in ("bus", "bus0", "bus1", "bus2", "bus3", "bus4"):
        if bcol in selector and bcol in base.columns:
            val = selector[bcol]
            if isinstance(val, str):
                val = [val]
            mask &= base[bcol].isin(list(val))

    for k, v in selector.items():
        if k in ("names", "carrier", "bus", "bus0", "bus1", "bus2", "bus3", "bus4"):
            continue
        if k in base.columns:
            if isinstance(v, (list, tuple, set)):
                mask &= base[k].isin(list(v))
            else:
                mask &= base[k].eq(v)

    return pd.Index(base.index[mask]).unique().tolist()


def _apply_patch_static(
    df: pd.DataFrame,
    col: str,
    scenario: str,
    names: list[str],
    op: str,
    value: float,
) -> None:
    """Apply a scalar patch to a static component table."""
    if not isinstance(df.index, pd.MultiIndex):
        if op == "set":
            df.loc[names, col] = value
        elif op == "scale":
            df.loc[names, col] = df.loc[names, col] * value
        elif op == "add":
            df.loc[names, col] = df.loc[names, col] + value
        else:
            raise ValueError(f"Unsupported op: {op}")
        return

    idx = pd.MultiIndex.from_product([[scenario], names], names=["scenario", "name"])
    if op == "set":
        df.loc[idx, col] = value
    elif op == "scale":
        df.loc[idx, col] = df.loc[idx, col] * value
    elif op == "add":
        df.loc[idx, col] = df.loc[idx, col] + value
    else:
        raise ValueError(f"Unsupported op: {op}")


def _apply_patch_timeseries(
    ts: pd.DataFrame,
    scenario: str,
    names: list[str],
    op: str,
    value: Any,
) -> None:
    """
    Apply a patch to a time series DataFrame.

    For stochastic networks, columns are MultiIndex (scenario, name).
    value can be:
    - scalar
    - array-like with length == len(ts.index)
    - DataFrame with columns matching names
    """
    if isinstance(ts.columns, pd.MultiIndex):
        if "scenario" in ts.columns.names:
            scenarios_avail = set(ts.columns.get_level_values("scenario"))
        else:
            scenarios_avail = set(ts.columns.get_level_values(0))

        if scenario not in scenarios_avail:
            logger.warning(
                "Timeseries patch: scenario '%s' not found in ts.columns; skipping.",
                scenario,
            )
            return

        if "name" in ts.columns.names:
            avail_names = set(ts.columns.get_level_values("name"))
        else:
            avail_names = set(ts.columns.get_level_values(1))
    else:
        avail_names = set(ts.columns)

    names = [n for n in pd.Index(names).unique().tolist() if n in avail_names]
    if not names:
        logger.warning(
            "Timeseries patch matched no existing columns for scenario '%s'; skipping.",
            scenario,
        )
        return

    if not isinstance(ts.columns, pd.MultiIndex):
        cols = names
    else:
        cols = pd.MultiIndex.from_product([[scenario], names], names=["scenario", "name"])

    if np.isscalar(value):
        if op == "set":
            ts.loc[:, cols] = value
        elif op == "scale":
            ts.loc[:, cols] = ts.loc[:, cols] * value
        elif op == "add":
            ts.loc[:, cols] = ts.loc[:, cols] + value
        else:
            raise ValueError(f"Unsupported op: {op}")
        return

    if isinstance(value, pd.DataFrame):
        v = value.reindex(ts.index)
        if v.isnull().values.any():
            raise ValueError("Provided DataFrame value has NaNs after reindexing to snapshots.")
        if not set(names).issubset(set(v.columns)):
            missing = sorted(set(names) - set(v.columns))
            raise ValueError(f"Provided DataFrame value missing columns: {missing}")
        v = v[names]

        if isinstance(ts.columns, pd.MultiIndex):
            v.columns = cols

        if op == "set":
            ts.loc[:, cols] = v.values
        elif op == "scale":
            ts.loc[:, cols] = ts.loc[:, cols].values * v.values
        elif op == "add":
            ts.loc[:, cols] = ts.loc[:, cols].values + v.values
        else:
            raise ValueError(f"Unsupported op: {op}")
        return

    arr = np.asarray(value)
    if arr.ndim != 1 or len(arr) != len(ts.index):
        raise ValueError(
            f"Array-like value must be 1D and match snapshots length ({len(ts.index)}); got shape {arr.shape}"
        )
    if op == "set":
        ts.loc[:, cols] = arr[:, None]
    elif op == "scale":
        ts.loc[:, cols] = ts.loc[:, cols].values * arr[:, None]
    elif op == "add":
        ts.loc[:, cols] = ts.loc[:, cols].values + arr[:, None]
    else:
        raise ValueError(f"Unsupported op: {op}")


def _normalize_component_table_name(component: str) -> str:
    """
    Normalize a user-facing component name to the corresponding network table name.

    Supported examples:
    - Generator -> generators
    - generators -> generators
    - Link -> links
    - Load -> loads
    """
    mapping = {
        "bus": "buses",
        "buses": "buses",
        "carrier": "carriers",
        "carriers": "carriers",
        "generator": "generators",
        "generators": "generators",
        "load": "loads",
        "loads": "loads",
        "line": "lines",
        "lines": "lines",
        "link": "links",
        "links": "links",
        "store": "stores",
        "stores": "stores",
        "storageunit": "storage_units",
        "storageunits": "storage_units",
        "storage_unit": "storage_units",
        "storage_units": "storage_units",
        "transformer": "transformers",
        "transformers": "transformers",
        "shuntimpedance": "shunt_impedances",
        "shuntimpedances": "shunt_impedances",
        "shunt_impedance": "shunt_impedances",
        "shunt_impedances": "shunt_impedances",
    }

    key = str(component).strip().lower().replace(" ", "").replace("-", "").replace(".", "")
    if key not in mapping:
        raise ValueError(f"Unsupported component '{component}'.")
    return mapping[key]


def _normalize_structured_scenario_spec(spec: Any) -> tuple[str | None, dict]:
    """
    Normalize a structured scenario specification.

    Accepted forms:
    - None
    - "scenario_name"
    - {"name": "scenario_name", "params": {...}}
    """
    if spec is None:
        return None, {}

    if isinstance(spec, str):
        return spec, {}

    if isinstance(spec, dict):
        name = spec.get("name")
        params = spec.get("params", {})
        if name is None:
            raise ValueError(
                "Structured scenario dict must contain key 'name'."
            )
        if not isinstance(params, dict):
            raise TypeError(
                f"structured_scenario.params must be a dict; got {type(params).__name__}"
            )
        return name, params

    raise TypeError(
        f"structured_scenario must be None, str, or dict; got {type(spec).__name__}"
    )


def _get_component_attr_tables(
    n: pypsa.Network,
    component: str,
    attribute: str,
) -> tuple[str, pd.DataFrame, pd.DataFrame | None]:
    """
    Return static and time-series tables for a component attribute.

    Returns
    -------
    table_name : str
        Base component table name, e.g. 'generators'
    static_df : pd.DataFrame
        Static component table, e.g. n.generators
    ts_df : pd.DataFrame | None
        Time-series attribute table, e.g. n.generators_t.p_max_pu, if it exists
    """
    table_name = _normalize_component_table_name(component)
    static_df = getattr(n, table_name)

    ts_df = None
    ts_container_name = f"{table_name}_t"
    if hasattr(n, ts_container_name):
        ts_container = getattr(n, ts_container_name)
        if hasattr(ts_container, attribute):
            ts_df = getattr(ts_container, attribute)

    return table_name, static_df, ts_df


def _available_timeseries_names(
    ts: pd.DataFrame,
    scenario: str | None = None,
) -> set[str]:
    """
    Return component names available in a time-series table.

    If scenario is provided and ts is stochastic, names are filtered to that scenario.
    """
    if not isinstance(ts.columns, pd.MultiIndex):
        return set(ts.columns.astype(str))

    if scenario is not None:
        if "scenario" in ts.columns.names:
            scenarios = ts.columns.get_level_values("scenario")
        else:
            scenarios = ts.columns.get_level_values(0)

        if "name" in ts.columns.names:
            names = ts.columns.get_level_values("name")
        else:
            names = ts.columns.get_level_values(1)

        return set(pd.Index(names[scenarios == scenario]).astype(str))

    if "name" in ts.columns.names:
        return set(ts.columns.get_level_values("name").astype(str))
    return set(ts.columns.get_level_values(-1).astype(str))


def _split_names_by_target(
    names: list[str],
    ts: pd.DataFrame | None,
    scenario: str | None = None,
) -> tuple[list[str], list[str]]:
    """
    Split matched component names into time-series-backed and static-only names.
    """
    if ts is None:
        return [], list(names)

    ts_names_avail = _available_timeseries_names(ts, scenario=scenario)
    names_ts = [name for name in names if str(name) in ts_names_avail]
    names_static = [name for name in names if str(name) not in ts_names_avail]
    return names_ts, names_static


def _validate_modify_rule(rule: Mapping[str, Any]) -> None:
    """Validate one generic modification rule."""
    required = {"component", "attribute", "operation", "value"}
    missing = sorted(required - set(rule))
    if missing:
        raise ValueError(f"Missing required keys in rule: {missing}")

    op = str(rule["operation"]).strip().lower()
    if op not in {"set", "scale", "add"}:
        raise ValueError(f"Unsupported operation '{op}'. Allowed: set, scale, add")

    target = str(rule.get("target", "auto")).strip().lower()
    if target not in {"auto", "static", "timeseries"}:
        raise ValueError(
            f"Unsupported target '{target}'. Allowed: auto, static, timeseries"
        )


def _apply_modify_components_rule(
    n: pypsa.Network,
    rule: Mapping[str, Any],
    scenario: str | None = None,
) -> None:
    """
    Apply one generic component modification rule.

    Rule format
    -----------
    {
      "component": "Generator",
      "attribute": "marginal_cost",
      "target": "auto" | "static" | "timeseries",
      "carrier": ["OCGT", "CCGT"],
      "operation": "scale" | "set" | "add",
      "value": 1.15
    }
    """
    rule = _ensure_dict(rule, "rule")
    _validate_modify_rule(rule)

    component = rule["component"]
    attribute = str(rule["attribute"])
    operation = str(rule["operation"]).strip().lower()
    target = str(rule.get("target", "auto")).strip().lower()
    value = rule["value"]

    table_name, static_df, ts_df = _get_component_attr_tables(
        n=n,
        component=component,
        attribute=attribute,
    )

    selector = {
        k: v
        for k, v in rule.items()
        if k not in {"component", "attribute", "target", "operation", "value"}
    }

    names = _select_names_from_component(n, table_name, selector)
    if not names:
        logger.warning(
            "modify_components rule matched no components. component=%s attribute=%s selector=%s",
            component,
            attribute,
            selector,
        )
        return

    if target == "timeseries" and ts_df is None:
        raise ValueError(
            f"Rule requested timeseries target for {component}.{attribute}, "
            f"but no time-series table exists."
        )

    if attribute == "p_max_pu" and target in {"timeseries", "auto"} and ts_df is not None:
        if operation != "scale":
            # Only forbid when the rule actually hits the timeseries branch
            names_ts, _ = _split_names_by_target(names, ts_df, scenario=scenario)
            if names_ts:
                raise ValueError(
                    "For timeseries p_max_pu only 'scale' is allowed."
                )

    if target == "static":
        if attribute not in static_df.columns:
            raise KeyError(f"Column not found in static table: {table_name}.{attribute}")
        _apply_patch_static(static_df, attribute, scenario, names, operation, value)
        logger.info(
            "Applied static modify_components rule on %s.%s to %s component(s)%s.",
            table_name,
            attribute,
            len(names),
            f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
        )
        return

    if target == "timeseries":
        _apply_patch_timeseries(ts_df, scenario, names, operation, value)
        logger.info(
            "Applied timeseries modify_components rule on %s_t.%s to %s component(s)%s.",
            table_name,
            attribute,
            len(names),
            f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
        )
        return

    # target == "auto"
    if ts_df is None:
        if attribute not in static_df.columns:
            raise KeyError(f"Column not found in static table: {table_name}.{attribute}")
        _apply_patch_static(static_df, attribute, scenario, names, operation, value)
        logger.info(
            "Applied auto->static modify_components rule on %s.%s to %s component(s)%s.",
            table_name,
            attribute,
            len(names),
            f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
        )
        return

    names_ts, names_static = _split_names_by_target(names, ts_df, scenario=scenario)

    if names_ts:
        _apply_patch_timeseries(ts_df, scenario, names_ts, operation, value)

    if names_static:
        if attribute not in static_df.columns:
            raise KeyError(f"Column not found in static table: {table_name}.{attribute}")
        _apply_patch_static(static_df, attribute, scenario, names_static, operation, value)

    logger.info(
        "Applied auto modify_components rule on %s.%s: %s timeseries-backed + %s static-only component(s)%s.",
        table_name,
        attribute,
        len(names_ts),
        len(names_static),
        f" for stochastic scenario '{scenario}'" if scenario is not None else " in deterministic mode",
    )


def _scenario_modify_components(
    n: pypsa.Network,
    scenario: str | None = None,
    config: dict | None = None,
) -> None:
    """
    Generic structured scenario applying one or more component modification rules.

    Expected config format
    ----------------------
    {
      "rules": [
        {
          "component": "Generator",
          "attribute": "marginal_cost",
          "target": "auto",
          "carrier": ["OCGT", "CCGT"],
          "operation": "scale",
          "value": 1.15
        }
      ]
    }
    """
    cfg = _ensure_dict(config, "config") if config is not None else {}
    rules = cfg.get("rules", None)
    if not isinstance(rules, list) or not rules:
        raise ValueError(
            "modify_components requires config['rules'] as a non-empty list."
        )

    for i, rule in enumerate(rules, start=1):
        logger.info(
            "Applying modify_components rule %s/%s%s.",
            i,
            len(rules),
            f" for stochastic scenario '{scenario}'" if scenario is not None else "",
        )
        _apply_modify_components_rule(n=n, rule=rule, scenario=scenario)

def _apply_patch_config_if_present(
    n: pypsa.Network,
    stochastic_param: dict,
) -> None:
    """
    Apply legacy patch-based stochastic configuration if patches are present.

    This keeps backward compatibility with the previous patch-based interface.
    Patches are only applied in stochastic mode because their values are keyed by scenario.
    """
    stoch = _merge_stochastic_param(stochastic_param)
    enabled = bool(stoch.get("enable", stoch.get("enabled", False)))
    if not enabled:
        return

    patches = stoch.get("patches", [])
    if not patches:
        logger.info("No legacy stochastic patches provided.")
        return

    def base_comp_from_table(table: str) -> str:
        return table[:-2] if table.endswith("_t") else table

    for i, patch in enumerate(patches, start=1):
        patch = _ensure_dict(patch, f"patch[{i}]")
        target = patch.get("target")
        if not isinstance(target, str) or "." not in target:
            raise ValueError(
                f"patch[{i}].target must be like 'generators.marginal_cost' or 'loads_t.p_set'"
            )

        table, attr = target.split(".", 1)
        selector = _ensure_dict(patch.get("selector", {}), f"patch[{i}].selector")
        op = patch.get("op", "set")
        values = _ensure_dict(patch.get("values", {}), f"patch[{i}].values")

        comp_for_selection = base_comp_from_table(table)
        names = _select_names_from_component(n, comp_for_selection, selector)
        if not names:
            logger.warning("patch[%s] matched no components; skipping. target=%s", i, target)
            continue

        logger.info("Applying patch[%s] target=%s op=%s matched=%s", i, target, op, len(names))

        if table.endswith("_t"):
            ts_container = getattr(n, table)
            ts = getattr(ts_container, attr)

            for sc, v in values.items():
                if isinstance(v, str) and v.endswith((".csv", ".parquet", ".pq")):
                    vp = Path(v)
                    if not vp.exists():
                        raise FileNotFoundError(f"patch[{i}] value file not found: {vp}")
                    if vp.suffix == ".csv":
                        dfv = pd.read_csv(vp, index_col=0, parse_dates=True)
                    else:
                        dfv = pd.read_parquet(vp)
                    _apply_patch_timeseries(ts, sc, names, op, dfv)
                else:
                    _apply_patch_timeseries(ts, sc, names, op, v)

        else:
            comp_df = getattr(n, table)
            if attr not in comp_df.columns:
                raise KeyError(f"patch[{i}] column not found: {table}.{attr}")

            for sc, v in values.items():
                if not np.isscalar(v):
                    raise ValueError(
                        f"patch[{i}] static patch values must be scalar; got {type(v).__name__}"
                    )
                _apply_patch_static(comp_df, attr, sc, names, op, float(v))


def apply_stochastic_config(
    n: pypsa.Network,
    config: dict,
    stochastic_param: dict,
    wildcards: Mapping[str, Any] | None = None,
) -> None:
    """
    Apply structured stochastic/deterministic scenarios and optional legacy patches.

    Behavior
    --------
    - If stochastic mode is enabled:
      * read scenario names from stochastic config
      * call n.set_scenarios(...)
      * dispatch one structured builder per scenario name
      * optionally apply legacy patch-based modifications

    - If stochastic mode is disabled:
      * detect a single structured scenario name from config or run context
      * apply the corresponding builder in deterministic mode
    """
    is_stochastic, active_names = _apply_structured_scenarios(
        n=n,
        config=config,
        stochastic_param=stochastic_param,
        wildcards=wildcards,
    )

    if is_stochastic:
        _apply_patch_config_if_present(n=n, stochastic_param=stochastic_param)

    if active_names:
        logger.info("Applied structured scenario(s): %s", active_names)
    else:
        logger.info("No structured scenario was applied.")


STRUCTURED_SCENARIOS = {
    "modify_components": _scenario_modify_components,
    "agriculture_full_electric": _scenario_agriculture_full_electric,
    "agriculture_machinery_full_oil": _scenario_agriculture_machinery_full_oil,
    "shipping_full_methanol": _scenario_shipping_full_methanol,
    "urban_heat_full_central": _scenario_urban_heat_full_central,
    "land_transport_linear_ev": _scenario_land_transport_linear_ev,
    "electricity_optimistic": _scenario_electricity_optimistic,
    "industry_h2": _scenario_industry_h2,
    "base": _scenario_base,
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_stochastic_network",
            opts="",
            clusters="adm",
            configfiles="config/test_stochastic_scenarios/config.yaml",
            sector_opts="",
            planning_horizons="2050",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    solve_opts = snakemake.params.solving["options"]
    np.random.seed(solve_opts.get("seed", 123))

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
        rolling_horizon=solve_opts.get("rolling_horizon", False),
    )

    apply_stochastic_config(
        n,
        config=snakemake.config,
        stochastic_param=snakemake.params.get("stochastic_scenarios", {}),
        wildcards=dict(snakemake.wildcards),
    )

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    with open(snakemake.output.config, "w", encoding="utf-8") as f:
        yaml.dump(
            n.meta,
            f,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )

    logger.info("Exported stochastic pre-solve network to %s", snakemake.output.network)