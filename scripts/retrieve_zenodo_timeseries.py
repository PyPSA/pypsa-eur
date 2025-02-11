# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract times series of wind power generation, solar power generation, hydropower inflow,
heating demand, and cooling demand that use both historical and projected climate variables.
Historical climate data are from the ERA5 reanalysis (1940-2023), while projected climate variables
are from three climate models of the CMIP5 EURO-CORDEX (2006-2100) and three representative concentration pathways (RCP 2.6, RCP 4.5 and RCP 8.5).
The time series are at country level for the EU27 countries, Balkan countries, the United Kingdom, Norway, and Switzerland.

Paper: Antonini et al., Weather- and climate-driven power supply and demand time series for power and energy system analyses https://doi.org/10.1038/s41597-024-04129-8
Zenodo repository: https://doi.org/10.5281/zenodo.13938926


"""

import logging
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_zenodo_timeseries")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    rcp = snakemake.params.rcp
    global_regional_models = snakemake.params.global_regional_models

    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    capacity_factor_urls = {
    "solar_capacity_factor": "https://zenodo.org/records/13938926/files/solar_capacity_factor.zip?download=1",
    "wind_capacity_factor": "https://zenodo.org/records/13938926/files/wind_capacity_factor.zip?download=1",
    "run_of_river_hydropower_inflow": "https://zenodo.org/records/13938926/files/run_of_river_hydropower_inflow.zip?download=1",
    "conventional_and_pumped_storage_hydropower_inflow": "https://zenodo.org/records/13938926/files/conventional_and_pumped_storage_hydropower_inflow.zip?download=1",
    }

    # Loop through the dictionary and download/extract each file
    for key, url in capacity_factor_urls.items():
        tarball_fn = Path(f"{rootpath}/data/zenodo_timeseries/{key}.zip")
        
        to_fn = Path(rootpath) / Path(snakemake.output[0])
        #to_fn = Path(f"{rootpath}/data/zenodo_timeseries/{key}/")

        logger.info(f"Downloading {key} data from '{url}'.")
        tarball_fn.parent.mkdir(parents=True, exist_ok=True)
        progress_retrieve(url, tarball_fn, disable=disable_progress) 

        logger.info(f"Extracting {key} capacity factor data.")
        with zipfile.ZipFile(tarball_fn, "r") as zip_ref:
            zip_ref.extractall(to_fn)

        logger.info(f"{key} data available in '{to_fn}'.")

        
        # Remove files that do not contain both rcp and global_regional_models in their names
        # ADB to do this does not work
        extracted_dir = to_fn / key
        for file in extracted_dir.iterdir():
            modified_name = file.name.replace('_', '-')
            if str(rcp) not in modified_name or global_regional_models not in modified_name:
                #logger.info(f"Removing file: {file}")
                file.unlink()
        
        # Remove the tarball (zip file)
        if tarball_fn.exists():
            logger.info(f"Removing zip file: {tarball_fn}")
            tarball_fn.unlink()

