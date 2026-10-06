# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve daily seawater temperature data from the Copernicus Marine Service.

Potential sea temperature (thetao) from the Global Ocean Physics Reanalysis is
downloaded for one year over European waters (12W to 42E, 33N to 72N) at 0.083
degree resolution and 5 to 15 m depth. The data provides the source
temperature of seawater heat pumps in district heating. For the test cutout, a
prepared subset is downloaded instead. The download requires Copernicus Marine
credentials configured for the copernicusmarine package.
"""

import logging
import os

import copernicusmarine
import requests

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "retrieve_seawater_temperature",
            year=2050,
        )

    # Configure logging and scenario
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    if snakemake.params.default_cutout == "be-03-2013-era5":
        logger.info("Retrieving test-cutout seawater temperature data.")

        url = snakemake.params.test_data_url

        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(snakemake.output.seawater_temperature, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        logger.info(
            f"Successfully downloaded test-cutout seawater temperature data to {snakemake.output.seawater_temperature}"
        )
    else:
        # Download seawater temperature data from Copernicus Marine Service
        # Dataset: Global Ocean Physics Reanalysis (daily, 0.083° resolution)
        # Variable: thetao (potential temperature in °C)
        # Spatial coverage: European waters (-12°W to 42°E, 33°N to 72°N)
        # Depth range: 5-15m (suitable for heat pump intake depths)
        logger.info(
            f"Downloading seawater temperature data for year {snakemake.wildcards.year}"
        )

        _ = copernicusmarine.subset(
            dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",  # Global ocean physics reanalysis
            start_datetime=f"{snakemake.wildcards.year}-01-01",
            end_datetime=f"{int(snakemake.wildcards.year)}-12-31",
            minimum_longitude=-12,  # Western European boundary
            maximum_longitude=42,  # Eastern European boundary
            minimum_latitude=33,  # Southern European boundary
            maximum_latitude=72,  # Northern European boundary
            variables=["thetao"],  # Potential temperature [°C]
            minimum_depth=5,  # Near-surface depth for heat pumps [m]
            maximum_depth=15,  # Near-surface depth for heat pumps [m]
            output_filename=snakemake.output.seawater_temperature,
        )

        # Verify successful download
        if not os.path.exists(snakemake.output.seawater_temperature):
            raise FileNotFoundError(
                f"Failed to retrieve seawater temperature data and save to {snakemake.output.seawater_temperature}. "
                f"One reason might be missing Copernicus Marine login info. "
                f"See the copernicusmarine package documentation for details."
            )

        logger.info(
            f"Successfully downloaded seawater temperature data to {snakemake.output.seawater_temperature}"
        )
