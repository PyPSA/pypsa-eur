# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Calculates dynamic line rating time series from base network.

Description
-----------

The rule :mod:`build_line_rating` calculates the line rating for transmission lines.
The line rating provides the maximal capacity of a transmission line considering the heat exchange with the environment.

The following heat gains and losses are considered:

- heat gain through resistive losses
- heat gain through solar radiation
- heat loss through radiation of the transmission line
- heat loss through forced convection with wind
- heat loss through natural convection


With a heat balance considering the maximum temperature threshold of the transmission line,
the maximal possible capacity factor "s_max_pu" for each transmission line at each time step is calculated.
"""

import logging
import re

import atlite
import geopandas as gpd
import numpy as np
import pypsa
import xarray as xr
from _helpers import configure_logging, get_snapshots, set_scenario_config
from dask.distributed import Client
from shapely.geometry import LineString as Line
from shapely.geometry import Point

logger = logging.getLogger(__name__)


def calculate_resistance(T, R_ref, T_ref: float | int = 293, alpha: float = 0.00403):
    """
    Calculates the resistance at other temperatures than the reference
    temperature.

    Parameters
    ----------
    T : Temperature at which resistance is calculated in [°C] or [K]
    R_ref : Resistance at reference temperature in [Ohm] or [Ohm/Per Length Unit]
    T_ref : Reference temperature in [°C] or [K]
    alpha: Temperature coefficient in [1/K]
        Defaults are:
            * T_ref : 20 °C
            * alpha : 0.00403 1/K

    Returns
    -------
    Resistance of at given temperature.
    """
    return R_ref * (1 + alpha * (T - T_ref))


def calculate_line_rating(
    n: pypsa.Network,
    cutout: atlite.Cutout,
    show_progress: bool = True,
    dask_kwargs: dict = None,
) -> xr.DataArray:
    """
    Calculates the maximal allowed power flow in each line for each time step
    considering the maximal temperature.

    Parameters
    ----------
    n : pypsa.Network object containing information on grid

    Returns
    -------
    xarray DataArray object with maximal power.
    """
    if dask_kwargs is None:
        dask_kwargs = {}

    logger.info("Calculating dynamic line rating.")
    relevant_lines = n.lines[~n.lines["underground"]].copy()
    buses = relevant_lines[["bus0", "bus1"]].values
    x = n.buses.x
    y = n.buses.y
    shapes = [Line([Point(x[b0], y[b0]), Point(x[b1], y[b1])]) for (b0, b1) in buses]
    shapes = gpd.GeoSeries(shapes, index=relevant_lines.index)
    if relevant_lines.r_pu.eq(0).all():
        # Overwrite standard line resistance with line resistance obtained from line type
        r_per_length = n.line_types["r_per_length"]
        R = (
            relevant_lines.join(r_per_length, on=["type"])["r_per_length"] / 1000
        )  # in meters
        # If line type with bundles is given retrieve number of conductors per bundle
        relevant_lines["n_bundle"] = (
            relevant_lines["type"]
            .where(relevant_lines["type"].str.contains("bundle"))
            .dropna()
            .apply(lambda x: int(re.findall(r"(\d+)-bundle", x)[0]))
        )
        # Set default number of bundles per line
        relevant_lines["n_bundle"] = relevant_lines["n_bundle"].fillna(1)
        R *= relevant_lines["n_bundle"]
        R = calculate_resistance(T=353, R_ref=R)
    Imax = cutout.line_rating(
        shapes,
        R,
        D=0.0218,
        Ts=353,
        epsilon=0.8,
        alpha=0.8,
        show_progress=show_progress,
        dask_kwargs=dask_kwargs,
    )
    line_factor = relevant_lines.eval("v_nom * n_bundle * num_parallel") / 1e3  # in mW
    return xr.DataArray(
        data=np.sqrt(3) * Imax * line_factor.values.reshape(-1, 1),
        attrs=dict(
            description="Maximal possible power in MW for given line considering line rating"
        ),
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_line_rating")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    show_progress = not snakemake.config["run"].get("disable_progressbar", True)
    show_progress = show_progress and snakemake.config["atlite"]["show_progress"]
    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
    else:
        client = None
    dask_kwargs = {"scheduler": client}

    n = pypsa.Network(snakemake.input.base_network)
    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=time)

    da = calculate_line_rating(n, cutout, show_progress, dask_kwargs)
    da.to_netcdf(snakemake.output[0])
