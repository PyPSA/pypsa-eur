# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Reproject NUTS2-level 1st-generation (1G) biofuel and perennial crop yields
(from ``build_perennials_yields_eurostat_average.py``) onto the clustered network
regions, for use by ``add_perennials()`` in ``prepare_sector_network.py``.

Yields are mapped with an area-weighted average over the NUTS2/region overlay.
NUTS2 regions without Eurostat data (AL, RS, BA, XK) take the yields of the
nearest NUTS2 region with data. Outputs one column per biomass class
(MWh of crop/ha/y for the 1G biofuel crops, t/ha/y for ``perennials``), indexed by
clustered region name.
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_biomass_potentials import (
    build_nuts2_shapes,
    convert_nuts2_to_regions,
    resolve_biomass_classes,
)

logger = logging.getLogger(__name__)


def impute_missing_values(df_nuts2, missing_shapes, yield_cols):
    """
    Fill missing yield values for NUTS2 regions not covered by Eurostat data
    (non-EU countries, small islands, city-states) by copying the values from
    the nearest NUTS2 region that does have valid (non-NaN) yields.

    Nearest is determined by centroid distance in an equal-area CRS (EPSG:3035).

    Parameters
    ----------
    df_nuts2 : gpd.GeoDataFrame
        NUTS2 geometries joined with yield columns; rows for regions without
        Eurostat coverage are entirely NaN in ``yield_cols``.
    missing_shapes : gpd.GeoDataFrame
        Geometries of the NUTS2 (or country-level substitute) regions to
        impute, indexed the same way as ``df_nuts2``.
    yield_cols : list of str
        Columns in ``df_nuts2`` to impute.

    Returns
    -------
    gpd.GeoDataFrame
        One row per entry in ``missing_shapes``, with ``yield_cols`` filled
        from the nearest valid NUTS2 neighbour.
    """
    df_valid = df_nuts2.dropna(subset=yield_cols).copy()

    valid_proj = df_valid.to_crs(3035)
    missing_proj = missing_shapes.to_crs(3035)

    valid_centroids = valid_proj.centroid
    missing_centroids = missing_proj.centroid

    imputed_rows = []

    for missing_id, c_geom in missing_centroids.items():
        dists = valid_centroids.distance(c_geom)
        nearest_nuts2 = dists.idxmin()

        attrs = df_valid.loc[nearest_nuts2, yield_cols]

        new_row = missing_shapes.copy().loc[[missing_id]]
        for col in yield_cols:
            new_row[col] = attrs[col]

        imputed_rows.append(new_row)

    return pd.concat(imputed_rows)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_perennials_yields",
            clusters="39",
            planning_horizons=2050,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    regions = gpd.read_file(snakemake.input.regions_onshore)
    nuts2 = build_nuts2_shapes(snakemake.input.nuts2, snakemake.input.country_shapes)

    yields = pd.read_csv(snakemake.input.perennials_yields_1G_biofuels, index_col=0)

    df_nuts2 = gpd.GeoDataFrame(nuts2.geometry).join(yields)

    # countries without Eurostat NUTS2 data take the yields of the nearest NUTS2 region
    missing_countries = df_nuts2.index.intersection(["AL", "RS", "BA", "XK"])
    if not missing_countries.empty:
        missing_shapes = df_nuts2.loc[missing_countries, ["geometry"]]
        yield_cols = [
            "Bioethanol barley, wheat, grain maize, oats, other cereals and rye",
            "Sugar from sugar beet",
            "Rape seed",
            "perennials",
        ]
        imputed_missing = impute_missing_values(df_nuts2, missing_shapes, yield_cols)
        df_nuts2 = pd.concat([df_nuts2.drop(index=missing_countries), imputed_missing])

    # convert nuts2 yields to regions
    df = convert_nuts2_to_regions(df_nuts2, regions, extensive=False)

    params = snakemake.params.biomass

    classes = resolve_biomass_classes(
        params["classes"], snakemake.config["sector"].get("perennials", False)
    )
    grouper = {v: k for k, vv in classes.items() for v in vv}

    # keep the column 'perennials' which is otherwise dropped
    for col in df.columns:
        if col not in grouper:
            grouper[col] = col

    df = df.T.groupby(grouper).sum().T

    df.index.name = "name"
    df.to_csv(snakemake.output.csv_file)
