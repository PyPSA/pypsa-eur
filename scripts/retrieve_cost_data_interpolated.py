# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve technology cost data from the technology-data repository.

The upstream dataset provides costs in 5-year increments. For intermediate
planning years, this script linearly interpolates between the nearest bracketing
anchor years and writes a synthetic costs_{planning_horizons}.csv.
"""

import logging
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)


def _anchor_years(year: int, *, min_year: int = 2020, max_year: int = 2050, step: int = 5) -> tuple[int, int]:
    """
    Return bracketing anchor years for the given target year.

    Parameters
    ----------
    year : int
        Target year.
    min_year, max_year : int
        Lower/upper bounds of the dataset coverage. Years outside the range are
        clamped to the nearest endpoint.
    step : int
        Spacing of anchor years (default: 5 years).

    Returns
    -------
    tuple[int, int]
        (before, after) anchor years. If year is on an anchor, returns (year, year).
    """
    y = int(year)

    if y <= min_year:
        return min_year, min_year
    if y >= max_year:
        return max_year, max_year
    if y % step == 0:
        return y, y

    before = (y // step) * step
    after = before + step
    return before, after


def _download(url: str, path: Path, *, timeout: int = 120) -> None:
    """
    Download a URL to a local file and fail fast on HTTP errors.
    """
    logger.info("Downloading %s", url)
    r = requests.get(url, stream=True, timeout=timeout)
    r.raise_for_status()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)


def _interpolate_costs_csv(
    path_before: Path,
    path_after: Path,
    year_before: int,
    year_after: int,
    year_target: int,
    path_out: Path,
) -> None:
    """
    Interpolate technology-data cost tables (long format) for a target year.

    The technology-data cost tables are expected to contain a numeric 'value'
    column and a set of identifier columns (e.g. technology/parameter/unit).
    The interpolation is performed on rows aligned by common identifier columns.

    Parameters
    ----------
    path_before, path_after : Path
        Input CSV paths for the two anchor years.
    year_before, year_after : int
        Anchor years.
    year_target : int
        Target year.
    path_out : Path
        Output CSV path.
    """
    y0 = int(year_before)
    y1 = int(year_after)
    yt = int(year_target)

    if y0 == y1:
        df = pd.read_csv(path_before)
        df.to_csv(path_out, index=False)
        return

    frac = (yt - y0) / (y1 - y0)

    df0 = pd.read_csv(path_before)
    df1 = pd.read_csv(path_after)

    if "value" not in df0.columns or "value" not in df1.columns:
        raise ValueError("Expected a 'value' column in technology-data cost tables.")

    # Prefer stable identifier columns if present; otherwise fall back to shared non-'value' columns.
    preferred_keys = ["technology", "parameter", "unit"]
    key_cols = [c for c in preferred_keys if c in df0.columns and c in df1.columns]
    if not key_cols:
        shared = [c for c in df0.columns if c in df1.columns and c != "value"]
        if not shared:
            raise ValueError("Could not determine common identifier columns to align cost tables.")
        key_cols = shared

    df0i = df0.set_index(key_cols)
    df1i = df1.set_index(key_cols)

    all_idx = df0i.index.union(df1i.index)

    v0 = pd.to_numeric(df0i.reindex(all_idx)["value"], errors="coerce")
    v1 = pd.to_numeric(df1i.reindex(all_idx)["value"], errors="coerce")

    # Interpolate where both endpoints exist; otherwise carry the available endpoint.
    v = v0 + frac * (v1 - v0)
    v = v.where(~(v0.isna() & ~v1.isna()), v1)
    v = v.where(~(v1.isna() & ~v0.isna()), v0)

    meta0 = df0i.reindex(all_idx).drop(columns=["value"], errors="ignore")
    meta1 = df1i.reindex(all_idx).drop(columns=["value"], errors="ignore")
    meta = meta0.combine_first(meta1)

    out = meta.copy()
    out["value"] = v.astype(float)

    path_out.parent.mkdir(parents=True, exist_ok=True)

    out_df = out.reset_index()

    # Reorder columns to match the original file column order (use df0 as template).
    template_cols = list(df0.columns)
    desired_cols = [c for c in template_cols if c in out_df.columns] + [
        c for c in out_df.columns if c not in template_cols
    ]
    out_df = out_df[desired_cols]

    out_df.to_csv(path_out, index=False)


if __name__ == "__main__":
    # Minimal logging setup compatible with Snakemake.
    log_file = None
    if "snakemake" in globals() and getattr(snakemake, "log", None):
        # snakemake.log may be a list-like container
        try:
            log_file = snakemake.log[0]
        except Exception:
            log_file = None

    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(levelname)s:%(name)s:%(message)s",
    )

    year = int(snakemake.wildcards.planning_horizons)
    base_url = str(snakemake.params.costs_url).rstrip("/")

    min_year = int(getattr(snakemake.params, "min_year", 2020))
    max_year = int(getattr(snakemake.params, "max_year", 2050))
    step = int(getattr(snakemake.params, "step", 5))

    before, after = _anchor_years(year, min_year=min_year, max_year=max_year, step=step)

    out_path = Path(snakemake.output.costs)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if before == after:
        url = f"{base_url}/costs_{before}.csv"
        _download(url, out_path)
        logger.info("Wrote %s (no interpolation needed).", out_path)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            f0 = tmp / f"costs_{before}.csv"
            f1 = tmp / f"costs_{after}.csv"

            _download(f"{base_url}/costs_{before}.csv", f0)
            _download(f"{base_url}/costs_{after}.csv", f1)

            _interpolate_costs_csv(f0, f1, before, after, year, out_path)

        logger.info("Wrote %s (interpolated between %s and %s).", out_path, before, after)

