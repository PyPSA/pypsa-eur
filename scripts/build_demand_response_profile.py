# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Builds per-device, per-horizon demand-response shiftable-load profiles
(MW, by bus) from the O'Reilly et al. (2022) openENTRANCE residential DR
dataset, already pre-processed by prepare_dr_data.py into NUTS2-level,
per-device, per-year CSVs (data/demand_response/achievable/).

For each device this redistributes its NUTS2-level "load" (= max possible
down-shift) and "max_increase" (= max possible up-shift) time series over
the model's buses, using population share within each NUTS2 region (from
nuts3_shapes.geojson) and area overlap of NUTS3 regions with bus regions
(redistribute_attribute, reused from build_electricity_demand_base.py).
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config
from scripts.build_electricity_demand_base import redistribute_attribute

logger = logging.getLogger(__name__)

# The O'Reilly dataset uses Eurostat NUTS codes, PyPSA-Eur config uses ISO2.
# Only Greece and the UK differ (mirrors prepare_dr_data.py).
ISO2_TO_EUROSTAT = {"GR": "EL", "GB": "UK"}

CET_TO_UTC = pd.Timedelta(hours=-1)


def read_device_variable(fn: str, eurostat_countries: list[str]) -> pd.DataFrame:
    """
    Read one prepare_dr_data.py output CSV: index (month, day, hour), one
    column per NUTS2 region, filtered to `eurostat_countries` (Eurostat-style
    prefixes, e.g. "EL" not "GR"). Column labels are whitespace-stripped.
    """
    df = pd.read_csv(fn, index_col=["month", "day", "hour"])
    df.columns = df.columns.str.strip()
    if eurostat_countries:
        df = df[[c for c in df.columns if c[:2] in eurostat_countries]]
    return df


def to_utc(df: pd.DataFrame, year: int, drop_leap_day: bool) -> pd.DataFrame:
    """
    Attach a real year to the generic (month, day, hour) index, treating it
    as CET (fixed +01:00, no DST), then shift to UTC. The one hour this
    leaves missing at the year boundary (CET Jan 1 00:00 of the *next* year)
    is filled by wrapping around to this same year's own Jan 1 00:00 CET row,
    since the data is a repeating representative year by construction.
    """
    month = df.index.get_level_values("month")
    day = df.index.get_level_values("day")
    hour = df.index.get_level_values("hour")
    if drop_leap_day:
        keep = ~((month == 2) & (day == 29))
        df, month, day, hour = df[keep], month[keep], day[keep], hour[keep]

    cet = pd.DatetimeIndex(
        pd.to_datetime({"year": year, "month": month, "day": day, "hour": hour})
    )
    out = df.set_axis(cet + CET_TO_UTC, axis=0).sort_index()
    wrap = out.iloc[[0]].set_axis([out.index[-1] + pd.Timedelta(hours=1)])
    return pd.concat([out, wrap])


def nuts2_bus_weights(
    nuts3_fn: str, onshore_regions_fn: str, nuts2_codes: list[str]
) -> pd.DataFrame:
    """
    (bus x NUTS2) weights redistributing each NUTS2 region over the model's
    buses: population share within the NUTS2 region (per NUTS3 sub-region),
    then area overlap of each NUTS3 region with the bus regions. Each column
    is rescaled to sum to exactly 1; a warning is logged only if the raw
    (pre-rescale) overlap covers less than 95% of the NUTS2 region's
    population (the rest falling outside the model's bus regions).
    """
    # nuts3_shapes.geojson has its own "name" column, which collides with
    # onshore_regions's "name" (bus name) during the overlay below and would
    # get silently suffixed away - rename it first (same fix as the existing
    # nuts3_distribution_keys() in build_electricity_demand_base.py).
    nuts3 = (
        gpd.read_file(nuts3_fn).rename(columns={"name": "nuts3_name"}).to_crs(epsg=3035)
    )
    regions = (
        gpd.read_file(onshore_regions_fn)
        .set_index("name", drop=False)
        .to_crs(epsg=3035)
    )

    weights = {}
    for code in nuts2_codes:
        group = nuts3[nuts3["level2"] == code].copy()
        if group.empty or group["pop"].sum() == 0:
            continue
        group["pop_share"] = group["pop"] / group["pop"].sum()
        redistributed = redistribute_attribute(
            group, regions.reset_index(drop=True), "pop_share"
        )
        total = redistributed.sum()
        if total == 0:
            continue
        if abs(total - 1) > 0.05:
            logger.warning(
                f"NUTS2 {code}: bus weights only cover {total:.1%} of the "
                "region before rescaling (rest falls outside the model's "
                "bus regions)"
            )
        weights[code] = redistributed / total  # rescale to sum to exactly 1

    return pd.DataFrame(weights).reindex(regions.index).fillna(0.0)


def report_region_mismatch(oreilly_codes: set[str], model_codes: set[str]) -> None:
    """
    Log NUTS2 codes that don't match between the O'Reilly data and the
    model's own NUTS3-derived regions (both already filtered to the model's
    countries).
    """
    missing_in_model = sorted(oreilly_codes - model_codes)
    missing_in_oreilly = sorted(model_codes - oreilly_codes)
    if missing_in_model:
        logger.warning(
            f"{len(missing_in_model)} NUTS2 codes in the O'Reilly data have "
            f"no match among the model's regions (dropped): {missing_in_model}"
        )
    if missing_in_oreilly:
        logger.warning(
            f"{len(missing_in_oreilly)} of the model's NUTS2 regions have no "
            f"O'Reilly data (left at zero): {missing_in_oreilly}"
        )


def build_demand_response_profile(
    device_csvs: dict[str, str],
    devices: dict,
    snapshots: dict,
    drop_leap_day: bool,
    nuts3_fn: str,
    onshore_regions_fn: str,
    countries: list[str],
) -> xr.Dataset:
    snap = get_snapshots(snapshots, drop_leap_day)
    year = int(snap[0].year)
    eurostat_countries = [ISO2_TO_EUROSTAT.get(c, c) for c in countries]

    raw = {}
    for dev in devices:
        raw[dev, "down"] = read_device_variable(
            device_csvs[f"load_{dev}"], eurostat_countries
        )
        raw[dev, "up"] = read_device_variable(
            device_csvs[f"max_increase_{dev}"], eurostat_countries
        )

    oreilly_codes = set().union(*(df.columns for df in raw.values()))

    nuts3 = gpd.read_file(nuts3_fn)
    model_codes = set(nuts3["level2"])
    if countries:
        model_codes = {c for c in model_codes if c[:2] in eurostat_countries}
    report_region_mismatch(oreilly_codes, model_codes)

    weights = nuts2_bus_weights(nuts3_fn, onshore_regions_fn, sorted(oreilly_codes))

    down_das, up_das = [], []
    for dev in devices:
        down_t = to_utc(raw[dev, "down"], year, drop_leap_day).reindex(
            snap, fill_value=0.0
        )
        up_t = to_utc(raw[dev, "up"], year, drop_leap_day).reindex(snap, fill_value=0.0)
        down_bus = down_t.reindex(columns=weights.columns, fill_value=0.0) @ weights.T
        up_bus = up_t.reindex(columns=weights.columns, fill_value=0.0) @ weights.T

        # weights.index inherits onshore_regions's index name ("name"); drop
        # it so xarray doesn't use it as the coordinate's dimension name.
        coords = {"time": snap, "bus": weights.index.rename(None)}
        down_das.append(
            xr.DataArray(down_bus, dims=["time", "bus"], coords=coords).assign_coords(
                device=dev
            )
        )
        up_das.append(
            xr.DataArray(up_bus, dims=["time", "bus"], coords=coords).assign_coords(
                device=dev
            )
        )

    return xr.Dataset(
        {
            "dr_max_down": xr.concat(down_das, dim="device"),
            "dr_max_up": xr.concat(up_das, dim="device"),
        }
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_demand_response_profile", horizon="2025")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params
    device_csvs = {
        f"{var}_{dev}": snakemake.input[f"{var}_{dev}"]
        for dev in params.devices
        for var in ["load", "max_increase"]
    }

    ds = build_demand_response_profile(
        device_csvs=device_csvs,
        devices=params.devices,
        snapshots=params.snapshots,
        drop_leap_day=params.drop_leap_day,
        nuts3_fn=snakemake.input.nuts3,
        onshore_regions_fn=snakemake.input.onshore_regions,
        countries=params.countries,
    )

    comp = dict(zlib=True, complevel=9, least_significant_digit=4)
    ds.to_netcdf(snakemake.output[0], encoding={v: comp for v in ds.data_vars})
