# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build hydrogen storage potentials of salt caverns per clustered region in TWh.

Potential cavern sites with their energy density from Caglayan et al. (2020)
are overlaid with the onshore and offshore regions. Each site's potential is
its energy density times its area, and the share of the site falling into a
region is credited to that region. Potentials are reported separately for
onshore sites, nearshore sites within 50 km of the coast and offshore sites.

!!! note "Data provenance"
    The site map was digitised from Figure 6 of the paper and scaled to the
    country totals of Figure 7, split by onshore, nearshore and offshore.

References
----------
- Caglayan et al. (2020), [Technical potential of salt caverns for hydrogen storage in Europe](https://doi.org/10.1016/j.ijhydene.2019.12.161)
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def concat_gdf(gdf_list, crs="EPSG:4326"):
    """
    Concatenate multiple geopandas dataframes with common coordinate reference
    system (crs).
    """
    return gpd.GeoDataFrame(pd.concat(gdf_list), crs=crs)


def load_bus_regions(onshore_path, offshore_path):
    """
    Load pypsa-eur on- and offshore regions and concat.
    """
    offshore_bus_regions = gpd.read_file(offshore_path)
    onshore_bus_regions = gpd.read_file(onshore_path)
    bus_regions = concat_gdf([offshore_bus_regions, onshore_bus_regions])
    bus_regions = bus_regions.dissolve(by="name", aggfunc="sum")

    return bus_regions


def area(gdf):
    """
    Returns area of GeoDataFrame geometries in square kilometers.
    """
    return gdf.to_crs(epsg=3035).area.div(1e6)


def salt_cavern_potential_by_region(caverns, regions):
    # calculate area of caverns shapes
    caverns["area_caverns"] = area(caverns)

    overlay = gpd.overlay(regions.reset_index(), caverns, keep_geom_type=True)

    # calculate share of cavern area inside region
    overlay["share"] = area(overlay) / overlay["area_caverns"]

    overlay["e_nom"] = overlay.eval(
        "capacity_per_area * share * area_caverns / 1000"
    )  # TWh

    return overlay.groupby(["name", "storage_type"]).e_nom.sum().unstack("storage_type")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_salt_cavern_potentials")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    fn_onshore = snakemake.input.onshore_regions
    fn_offshore = snakemake.input.offshore_regions

    regions = load_bus_regions(fn_onshore, fn_offshore)

    caverns = gpd.read_file(snakemake.input.salt_caverns)  # GWh/sqkm

    caverns_regions = salt_cavern_potential_by_region(caverns, regions)

    caverns_regions.to_csv(snakemake.output.h2_cavern_potential)
