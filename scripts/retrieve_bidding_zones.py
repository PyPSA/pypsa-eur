# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieves bidding zone shape files from two sources. `electricitymaps-contrib` provides shape files for all the zones on a global level. `entsoe-py` provides country level shape files which are concatenated into one file. The `electricitymaps-contrib` data is preferred, but the Italian bidding zones from `entsoe-py` are more accurate.

Outputs
-------

- ``data/busshapes/bidding_zones_electricitymaps.geojson``:
- ``data/busshapes/bidding_zones_entsoepy.geojson``:
"""

from urllib.error import HTTPError

import entsoe
import geopandas as gpd
import pandas as pd


def load_bidding_zones_from_entsoepy() -> gpd.GeoDataFrame:
    """
    Load bidding zone geometries from entsoe-py GeoJSON files with disk caching.

    Returns:
        GeoDataFrame: Contains geometries for all available bidding zones
    """
    # If not in cache or cache disabled, load from source
    print("Downloading bidding zones...")
    gdfs: list[gpd.GeoDataFrame] = []
    for area in entsoe.Area:
        name = area.name
        try:
            url = f"https://raw.githubusercontent.com/EnergieID/entsoe-py/c03c604af36ef92e8ef6ee89dc57c56ca5e1dbac/entsoe/geo/geojson/{name}.geojson"
            gdfs.append(gpd.read_file(url))
        except HTTPError:
            continue

    shapes = pd.concat(gdfs, ignore_index=True)  # type: ignore

    return shapes


def load_bidding_zones_from_electricitymaps() -> gpd.GeoDataFrame:
    """
    Load bidding zone geometries from electricitymaps-contrib repository.

    Returns:
        GeoDataFrame: Contains geometries for all available bidding zones
    """
    url = "https://raw.githubusercontent.com/electricitymaps/electricitymaps-contrib/v1.238.0/web/geo/world.geojson"
    return gpd.read_file(url)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_bidding_zones")

    bidding_zones = load_bidding_zones_from_entsoepy()
    bidding_zones.to_file(snakemake.output.file_entsoepy)

    bidding_zones = load_bidding_zones_from_electricitymaps()
    bidding_zones.to_file(snakemake.output.file_electricitymaps)
