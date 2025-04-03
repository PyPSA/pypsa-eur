"""
Retrieve seawater temperature from the Copernicus Marine Service within the cutout bounds and save it to a NetCDF file.

Raises
------
ValueError
    If min/max cutout bounds are outside Europe or min exceed max bounds.

Inputs
------
- `cutouts/{cutout}.nc`: The cutout file containing the bounds of the area to retrieve the seawater data for.

Outputs
-------
- `data/{run_name}/seawater_temperature.nc`: The seawater data within the cutout bounds. Contains
"""

import logging

import atlite
import copernicusmarine
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_brownfield",
            clusters="39",
            opts="",
            ll="vopt",
            sector_opts="",
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    cutout = atlite.Cutout(snakemake.input.cutout)

    minimum_longitude = cutout.bounds[0]
    minimum_latitude = cutout.bounds[1]
    maximum_longitude = cutout.bounds[2]
    maximum_latitude = cutout.bounds[3]

    # sanity-check bounds
    if (
        minimum_longitude > minimum_longitude
        or minimum_latitude > maximum_latitude
        or minimum_longitude < -25
        or maximum_longitude > 45
        or minimum_latitude < 30
        or maximum_latitude > 75
    ):
        raise ValueError("Invalid bounds")

    data = copernicusmarine.subset(
        dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
        start_datetime=snakemake.params.snapshots["start"],
        end_datetime=snakemake.params.snapshots["end"],
        minimum_longitude=cutout.bounds[0],
        maximum_longitude=cutout.bounds[2],
        minimum_latitude=cutout.bounds[1],
        maximum_latitude=cutout.bounds[3],
        minimum_depth=40,
        maximum_depth=50,
        variables=["thetao"],
        output_filename=snakemake.output[0],
        disable_progress_bar=True,
    )
