# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Add transmission projects and DLR to the network.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def attach_transmission_projects(
    n: pypsa.Network, transmission_projects: list[str]
) -> None:
    logger.info("Adding transmission projects to network.")
    for path in transmission_projects:
        path = Path(path)
        df = pd.read_csv(path, index_col=0, dtype={"bus0": str, "bus1": str})
        if df.empty:
            continue
        if "new_buses" in path.name:
            n.add("Bus", df.index, **df)
        elif "new_lines" in path.name:
            n.add("Line", df.index, **df)
        elif "new_links" in path.name:
            n.add("Link", df.index, **df)
        elif "adjust_lines" in path.name:
            n.lines.update(df)
        elif "adjust_links" in path.name:
            n.links.update(df)


def attach_line_rating(
    n: pypsa.Network,
    rating: pd.DataFrame,
    s_max_pu: float,
    correction_factor: float,
    max_voltage_difference: float | bool,
    max_line_rating: float | bool,
) -> None:
    logger.info("Attaching dynamic line rating to network.")
    # TODO: Only considers overhead lines
    n.lines_t.s_max_pu = (rating / n.lines.s_nom[rating.columns]) * correction_factor
    if max_voltage_difference:
        x_pu = (
            n.lines.type.map(n.line_types["x_per_length"])
            * n.lines.length
            / (n.lines.v_nom**2)
        )
        # need to clip here as cap values might be below 1
        # -> would mean the line cannot be operated at actual given pessimistic ampacity
        s_max_pu_cap = (
            np.deg2rad(max_voltage_difference) / (x_pu * n.lines.s_nom)
        ).clip(lower=1)
        n.lines_t.s_max_pu = n.lines_t.s_max_pu.clip(
            lower=1, upper=s_max_pu_cap, axis=1
        )
    if max_line_rating:
        n.lines_t.s_max_pu = n.lines_t.s_max_pu.clip(upper=max_line_rating)
    n.lines_t.s_max_pu *= s_max_pu


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_transmission_projects_and_dlr")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params

    n = pypsa.Network(snakemake.input.network)

    if params["transmission_projects"]["enable"]:
        attach_transmission_projects(n, snakemake.input.transmission_projects)

    if params["dlr"]["activate"]:
        rating = xr.open_dataarray(snakemake.input.dlr).to_pandas().transpose()

        s_max_pu = params["s_max_pu"]
        correction_factor = params["dlr"]["correction_factor"]
        max_voltage_difference = params["dlr"]["max_voltage_difference"]
        max_line_rating = params["dlr"]["max_line_rating"]

        attach_line_rating(
            n,
            rating,
            s_max_pu,
            correction_factor,
            max_voltage_difference,
            max_line_rating,
        )

    n.export_to_netcdf(snakemake.output[0])
