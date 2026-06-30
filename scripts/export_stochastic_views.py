# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Export deterministic views from a stochastic PyPSA network solution.

Two modes:
- average: create __avg.nc where scenario-dependent tables are probability-weighted averages
- scenario: create __sc-{scenario}.nc where scenario-dependent tables are sliced at a given scenario

This script is called by Snakemake rules:
- export_stochastic_average
- export_stochastic_scenario
"""

import logging
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import pypsa
import yaml

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]  # points to /dati/pampado/pypsa-eur
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def _find_scenario_level(mi: pd.MultiIndex, scenario_names: list[str]) -> int | None:
    """
    Return the level index in a MultiIndex that corresponds to scenarios.

    Works even if the level name is None, by matching the set of values against scenario_names.
    """
    if not isinstance(mi, pd.MultiIndex):
        return None

    scen_set = set(map(str, scenario_names))
    for lvl in range(mi.nlevels):
        vals = pd.Index(mi.get_level_values(lvl).unique()).map(str)
        if len(vals) > 0 and set(vals).issubset(scen_set):
            return lvl
    return None


def read_probabilities(path: str) -> dict[str, float]:
    """Read stochastic scenario probabilities from a YAML file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Stochastic scenarios YAML not found: {p}")

    with p.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    scenarios = data.get("scenarios", data) or {}
    if not isinstance(scenarios, dict) or not scenarios:
        raise ValueError("Stochastic scenario file must contain a non-empty mapping.")

    probs = {}

    for name, spec in scenarios.items():
        if np.isscalar(spec):
            probs[str(name)] = float(spec)
            continue

        if isinstance(spec, dict):
            if "prob" in spec:
                probs[str(name)] = float(spec["prob"])
            elif "probability" in spec:
                probs[str(name)] = float(spec["probability"])
            else:
                raise ValueError(
                    f"Scenario '{name}' is missing 'prob' or 'probability'."
                )
            continue

        raise TypeError(
            f"Invalid scenario specification for '{name}': "
            f"expected scalar or dict, got {type(spec).__name__}."
        )

    total = sum(probs.values())
    if total <= 0:
        raise ValueError("Scenario probabilities must sum to a positive value.")

    if not np.isclose(total, 1.0, atol=1e-9):
        logger.warning(
            "Scenario probabilities sum to %.12g; normalizing to 1.0.",
            total,
        )
        probs = {k: v / total for k, v in probs.items()}

    return probs


def _is_scenario_index(idx: pd.Index) -> bool:
    return isinstance(idx, pd.MultiIndex) and ("scenario" in idx.names) and ("name" in idx.names)


def _is_scenario_columns(cols: pd.Index) -> bool:
    return isinstance(cols, pd.MultiIndex) and ("scenario" in cols.names) and ("name" in cols.names)


def _slice_static(df: pd.DataFrame, scenario: str, scenario_names: list[str]) -> pd.DataFrame:
    if getattr(df, "empty", False):
        return df
    if not isinstance(df.index, pd.MultiIndex):
        return df

    lvl = None
    if "scenario" in df.index.names:
        lvl = df.index.names.index("scenario")
    else:
        lvl = _find_scenario_level(df.index, scenario_names)

    if lvl is None:
        return df

    available = pd.Index(df.index.get_level_values(lvl).unique()).map(str)

    if str(scenario) in set(available):
        out = df.xs(scenario, level=lvl, drop_level=True)
        return out

    # fallback: scenario-invariant
    if len(available) == 1:
        s0 = df.index.get_level_values(lvl).unique()[0]
        logger.warning("Static table has only scenario=%r; using it for scenario=%r.", s0, scenario)
        return df.xs(s0, level=lvl, drop_level=True)

    raise KeyError(f"Scenario {scenario!r} not found. Available: {available.tolist()}")





def _expected_static(df: pd.DataFrame, probs: Dict[str, float]) -> pd.DataFrame:
    """Probability-weighted expected value for static tables (rare, but safe)."""
    if not _is_scenario_index(df.index):
        return df

    scen = df.index.get_level_values("scenario")
    w = scen.map(probs).to_numpy()

    num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    other_cols = [c for c in df.columns if c not in num_cols]

    out_parts = []

    if other_cols:
        out_other = df[other_cols].groupby(level="name").first()
        out_parts.append(out_other)

    if num_cols:
        out_num = df[num_cols].mul(w, axis=0).groupby(level="name").sum()
        out_parts.append(out_num)

    out = pd.concat(out_parts, axis=1) if out_parts else df.groupby(level="name").first()
    out = out.reindex(columns=[c for c in df.columns if c in out.columns])
    out.index.name = "name"
    return out


def _slice_timeseries(ts: pd.DataFrame | pd.Series, scenario: str, scenario_names: list[str]):
    if getattr(ts, "empty", False):
        return ts

    # We expect DataFrame with snapshots index and columns as names (or MultiIndex).
    if isinstance(ts, pd.Series):
        # Rare but handle: MultiIndex on index
        if not isinstance(ts.index, pd.MultiIndex):
            return ts

        lvl = None
        if "scenario" in ts.index.names:
            lvl = ts.index.names.index("scenario")
        else:
            lvl = _find_scenario_level(ts.index, scenario_names)

        if lvl is None:
            return ts

        available = pd.Index(ts.index.get_level_values(lvl).unique()).map(str)
        if str(scenario) in set(available):
            return ts.xs(scenario, level=lvl, drop_level=True)

        if len(available) == 1:
            s0 = ts.index.get_level_values(lvl).unique()[0]
            logger.warning("Series has only scenario=%r; using it for scenario=%r.", s0, scenario)
            return ts.xs(s0, level=lvl, drop_level=True)

        raise KeyError(f"Scenario {scenario!r} not found in Series. Available: {available.tolist()}")

    # DataFrame case
    if not isinstance(ts.columns, pd.MultiIndex):
        return ts

    lvl = None
    if "scenario" in ts.columns.names:
        lvl = ts.columns.names.index("scenario")
    else:
        lvl = _find_scenario_level(ts.columns, scenario_names)

    if lvl is None:
        return ts

    available = pd.Index(ts.columns.get_level_values(lvl).unique()).map(str)
    if str(scenario) in set(available):
        out = ts.xs(scenario, axis=1, level=lvl, drop_level=True)
        return out

    if len(available) == 1:
        s0 = ts.columns.get_level_values(lvl).unique()[0]
        logger.warning("Time series has only scenario=%r; using it for scenario=%r.", s0, scenario)
        return ts.xs(s0, axis=1, level=lvl, drop_level=True)

    raise KeyError(f"Scenario {scenario!r} not found in columns. Available: {available.tolist()}")



def _expected_timeseries(ts: pd.DataFrame, probs: dict[str, float]) -> pd.DataFrame:
    """Probability-weighted expected value for time series tables."""
    # If it's empty, expected is empty; just make it deterministic-consistent
    if ts is None or not isinstance(ts, pd.DataFrame) or ts.shape[1] == 0:
        if isinstance(ts, pd.DataFrame) and isinstance(ts.columns, pd.MultiIndex):
            # Drop scenario level if present (even if empty) to avoid scenario columns in deterministic view
            col_names = list(ts.columns.names)
            if "scenario" in col_names:
                sc_level = "scenario"
            else:
                sc_level = col_names[0] if col_names and col_names[0] is not None else 0

            # If there is at least a second level, keep that as "name"; else just empty Index
            if isinstance(ts.columns, pd.MultiIndex) and ts.columns.nlevels >= 2:
                # Create an empty Index with proper name
                ts = ts.copy()
                ts.columns = pd.Index([], name="name")
            else:
                ts = ts.copy()
                ts.columns = pd.Index([], name="name")
        return ts

    if not _is_scenario_columns(ts.columns):
        return ts

    # --- Identify scenario level name robustly ---
    col_names = list(ts.columns.names)
    if "scenario" in col_names:
        sc_level = "scenario"
        name_level = "name" if "name" in col_names else None
    else:
        sc_level = col_names[0] if col_names and col_names[0] is not None else 0
        name_level = col_names[1] if len(col_names) > 1 and col_names[1] is not None else 1

    scenarios_in_table = pd.Index(ts.columns.get_level_values(sc_level)).unique().tolist()

    def _norm(x) -> str:
        return str(x).strip()

    available = {_norm(s): s for s in scenarios_in_table}

    # --- Normalize/unwrap probs ---
    p = probs or {}
    if isinstance(p, dict) and "scenarios" in p and isinstance(p["scenarios"], dict):
        p = p["scenarios"]

    weights: dict[object, float] = {}
    for k, v in (p.items() if isinstance(p, dict) else []):
        if isinstance(v, dict):
            if "probability" in v:
                v = v["probability"]
            elif "p" in v:
                v = v["p"]
            else:
                continue
        if v is None:
            continue

        nk = _norm(k)
        if nk not in available:
            continue

        try:
            w = float(v)
        except Exception:
            continue

        weights[available[nk]] = w

    # If the table has scenario-columns but none match, do NOT crash if table is effectively empty.
    # Here ts has columns, so it's meaningful: raise to catch config/name mismatches.
    if not weights:
        raise ValueError(
            "No matching scenarios found in time series columns for expected export. "
            f"Available scenarios in table: {scenarios_in_table}. "
            f"Provided probs keys: {list(p.keys()) if isinstance(p, dict) else type(p)}."
        )

    wsum = sum(weights.values())
    if wsum <= 0:
        raise ValueError(f"Scenario probabilities sum to {wsum}, cannot compute expected value.")
    weights = {sc: w / wsum for sc, w in weights.items()}

    out = None
    for sc, w in weights.items():
        part = ts.xs(sc, level=sc_level, axis=1, drop_level=True) * w
        out = part if out is None else out.add(part, fill_value=0.0)

    if isinstance(out.columns, pd.Index):
        out.columns.name = "name"

    return out



def _iter_static_tables(n: pypsa.Network) -> Tuple[str, pd.DataFrame]:
    """Yield (attr_name, df) for all DataFrame attributes that look like component tables."""
    for attr in dir(n):
        if attr.startswith("_"):
            continue
        try:
            obj = getattr(n, attr)
        except Exception:
            continue
        if isinstance(obj, pd.DataFrame):
            yield attr, obj


def _iter_timeseries_tables(n: pypsa.Network):
    """
    Yield (tname, field, obj) for ALL time-dependent tables present in <something>_t containers.

    Notes
    -----
    - This iterates directly over Network attributes ending with '_t' (e.g. buses_t, loads_t, ...)
      to avoid relying on n.components iteration behavior (which can yield objects, not strings).
    - It captures both input time series and outputs (e.g. marginal_price, p, ...).
    """
    for tname in dir(n):
        if not tname.endswith("_t"):
            continue
        if not hasattr(n, tname):
            continue

        container = getattr(n, tname)
        if container is None:
            continue

        # PyPSA *_t containers are usually dict-like (AttrDict / dict)
        if hasattr(container, "items"):
            for field, obj in container.items():
                if isinstance(obj, (pd.DataFrame, pd.Series)):
                    yield tname, field, obj
        else:
            # Fallback: attribute-based container
            for field in dir(container):
                if field.startswith("_"):
                    continue
                obj = getattr(container, field)
                if isinstance(obj, (pd.DataFrame, pd.Series)):
                    yield tname, field, obj




def build_view(n: pypsa.Network, mode: str, probs: Dict[str, float], scenario: str | None) -> pypsa.Network:
    """Return a deterministic-view network."""
    m = n.copy()

    if mode not in {"average", "scenario"}:
        raise ValueError(f"Unknown mode={mode}")

    if mode == "scenario" and not scenario:
        raise ValueError("Scenario mode requires a scenario name.")

    # Static tables
    for attr, df in _iter_static_tables(m):
        if mode == "scenario":
            new_df = _slice_static(df, scenario, m.scenarios.tolist())
        else:
            new_df = _expected_static(df, probs)
        try:
            setattr(m, attr, new_df)
        except Exception:
            # Some attributes might be read-only; ignore safely.
            pass

    # Time series tables
    for tname, field, df in _iter_timeseries_tables(m):
        if mode == "scenario":
            new_ts = _slice_timeseries(df, scenario, m.scenarios.tolist())
        else:
            new_ts = _expected_timeseries(df, probs)
        try:
            setattr(getattr(m, tname), field, new_ts)
        except Exception:
            pass

    return m


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("export_stochastic_average",
                                    configfiles=["config/test_stoch/config.yaml"],
                                    clusters='adm',
                                    opts='',
                                    sector_opts='',
                                    planning_horizons='2050',
                                    stoch_scenario='high')

    configure_logging(snakemake)
    n = pypsa.Network(snakemake.input.network)

    probs = read_probabilities(snakemake.params.scenarios_file)

    mode = getattr(snakemake.params, "mode", "average")
    scenario = getattr(snakemake.params, "scenario", None)

    logger.info(f"Export mode={mode} scenario={scenario}")
    m = build_view(n, mode=mode, probs=probs, scenario=scenario)

    # Preserve meta if present
    try:
        m.meta = getattr(n, "meta", None)
    except Exception:
        pass

    out = snakemake.output.network

    m.export_to_netcdf(out)
    logger.info(f"Exported deterministic view to {out}")
