# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create a weather data cutout with [atlite](https://atlite.readthedocs.io/en/latest/).

The cutout holds gridded hourly weather variables for the configured bounding
box and period, downloaded from the
[ERA5](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)
reanalysis and optionally amended with
[SARAH-3](https://wui.cmsaf.eu/safira/action/viewProduktSearch) satellite-based
radiation observations. Downloading ERA5 requires the `cdsapi` package and a
registered [Copernicus Climate Data Store API key](https://cds.climate.copernicus.eu/how-to-api).
It is the weather input of the renewable profile, hydro inflow and temperature
profile rules downstream. See the
[atlite documentation](https://atlite.readthedocs.io/en/latest/examples/create_cutout.html)
for details on creating cutouts.

| Field | Dimensions | Unit | Description |
| --- | --- | --- | --- |
| pressure | time, y, x | Pa | Surface pressure |
| temperature | time, y, x | K | Air temperature 2 m above the surface |
| soil temperature | time, y, x | K | Soil temperature between 1 m and 3 m depth (layer 4) |
| influx_toa | time, y, x | W/m2 | Top of atmosphere incident solar radiation |
| influx_direct | time, y, x | W/m2 | Total sky direct solar radiation at surface |
| influx_diffuse | time, y, x | W/m2 | Diffuse solar radiation at surface (downward minus direct) |
| albedo | time, y, x | - | Share of downward solar radiation reflected by the surface, between 0 and 1 |
| runoff | time, y, x | m | Surface runoff (volume per area) |
| wnd100m | time, y, x | m/s | Wind speed at 100 m |
| roughness | y, x | m | Forecast surface roughness length |
| height | y, x | m | Surface elevation above sea level |

![](../img/era5.png)

A SARAH-3 cutout amends the ERA5 fields `temperature`, `influx_toa`,
`influx_direct`, `influx_diffuse` and `albedo`.

![](../img/sarah.png)
"""

import logging

import atlite

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_cutout", cutout="europe-2013-sarah3-era5")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    cutout_params = snakemake.params.cutouts[snakemake.wildcards.cutout]
    cutout_params["time"] = slice(*cutout_params["time"])
    cutout_params["x"] = slice(*cutout_params["x"])
    cutout_params["y"] = slice(*cutout_params["y"])
    prepare_kwargs = cutout_params.pop("prepare_kwargs", {})

    logger.info(f"Creating cutout with parameters {cutout_params}.")
    cutout = atlite.Cutout(snakemake.output[0], **cutout_params)

    logger.info(f"Preparing cutout the cutout with parameters {prepare_kwargs}.")
    cutout.prepare(**prepare_kwargs)
