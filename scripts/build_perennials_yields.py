# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Reproject NUTS2-level 1st-generation (1G) biofuel and perennial crop YIELDS
(from ``build_perennials_yields_eurostat_average.py``) onto the clustered network
regions, for use by ``add_perennials()`` in ``prepare_sector_network.py``.

Output units, NOT an area or a potential
-----------------------------------------
The output CSV holds per-hectare yield RATES, not land area or a CO2
potential: MWh/ha/y for the 1G biofuel crop columns (cereals, sugar beet,
rapeseed) and t/ha/y for the ``perennials`` column. The land area available
for conversion to perennial grasses, and the resulting CO2 sequestration
potential, are only derived later - see "Downstream usage" below.

Method
------
NUTS2 yields are mapped to clustered regions via an area-weighted average
over the NUTS2/region overlay (NUTS2 geometries intersected with cluster
region geometries, weighted by intersection area) - see
``convert_nuts2_to_regions_yields()``. NUTS2 regions not covered by
Eurostat data (non-EU countries, small islands, city-states) are first
filled from the nearest valid NUTS2 centroid - see
``impute_missing_values()``. Resulting columns are then grouped into
biomass classes (``resolve_biomass_classes()``) to match the class
structure used by ``build_biomass_potentials.py``.

Downstream usage
-----------------
``add_perennials()`` divides the ENSPRESO 1G biomass potential (MWh/y, an
extensive quantity from ``build_biomass_potentials.py``) by this script's
1G yields (MWh/ha/y) to back out the land area (ha) available for
conversion to perennial grasses, then converts that area into a CO2
sequestration potential via ``perennials.sequestration_co2`` (tCO2e/ha/y).

Outputs a single CSV with one column per crop class (cereals, sugar beet,
rapeseed, perennials) indexed by clustered region name.
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import (
    configure_logging,
    resolve_biomass_classes,
    set_scenario_config,
)
from scripts.build_biomass_potentials import build_nuts2_shapes

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
    # Keep only rows that have valid yields (drop NaN rows!)
    df_valid = df_nuts2.dropna(subset=yield_cols).copy()

    # Project to 3035
    valid_proj = df_valid.to_crs(3035)
    missing_proj = missing_shapes.to_crs(3035)

    valid_centroids = valid_proj.centroid
    missing_centroids = missing_proj.centroid

    imputed_rows = []

    for missing_id, c_geom in missing_centroids.items():

        # Distance to only VALID NUTS2 rows
        dists = valid_centroids.distance(c_geom)
        nearest_nuts2 = dists.idxmin()

        # Copy yields
        attrs = df_valid.loc[nearest_nuts2, yield_cols]

        # Create new row
        new_row = missing_shapes.copy().loc[[missing_id]]
        for col in yield_cols:
            new_row[col] = attrs[col]

        imputed_rows.append(new_row)

    return pd.concat(imputed_rows)


def convert_nuts2_to_regions_yields(df_nuts2, regions, yield_cols=None):
    """
    Convert NUTS2-level yields (intensive, e.g. MWh/ha/y or t/ha/y) to
    PyPSA-Eur clustered regions via an area-weighted average over the
    NUTS2/region overlay:

        y_n = Σ_i (y_i * A_i∩n) / Σ_i A_i∩n

    Unlike an extensive quantity (a total, e.g. MWh), a yield is a rate and
    must be area-weighted-averaged rather than redistributed by area share
    - see ``convert_nuts2_to_regions`` in ``build_biomass_potentials.py`` for
    the extensive-quantity equivalent.

    Only NUTS2 rows with non-NaN yields are used. Regions with no
    overlapping valid NUTS2 get NaN.

    Parameters
    ----------
    df_nuts2 : gpd.GeoDataFrame
        NUTS2 geometries joined with yield columns.
    regions : gpd.GeoDataFrame
        PyPSA-Eur clustered onshore regions, with a ``name`` column.
    yield_cols : list of str, optional
        Columns in ``df_nuts2`` to convert. Defaults to all columns except
        ``geometry`` and ``NUTS_ID``.

    Returns
    -------
    pd.DataFrame
        Area-weighted-average yields indexed by region name, one row per
        entry in ``regions``.
    """

    nuts = df_nuts2.copy()
    regs = regions.copy()

    # Ensure ID column for NUTS2
    nuts["NUTS_ID"] = nuts.index

    # Identify yield columns if not given
    if yield_cols is None:
        yield_cols = nuts.columns.difference(["geometry", "NUTS_ID"])

    # 1) Keep only NUTS2 rows that have at least one non-NaN yield
    nuts_valid = nuts.dropna(subset=yield_cols, how="all")

    # 2) Reproject to equal-area CRS for areas
    nuts_valid = nuts_valid.to_crs(3035)
    regs = regs.to_crs(3035)

    # 3) Overlay: intersection between regions and valid NUTS2
    overlay = gpd.overlay(regs, nuts_valid, keep_geom_type=True)

    # Return empty if nothing overlaps
    if overlay.empty:
        return pd.DataFrame(index=regs["name"].values, columns=yield_cols, dtype=float)

    # 4) Area of intersections (m² → km²)
    overlay["area_intersection"] = overlay.geometry.area

    # Optional: drop absurdly tiny slivers
    overlay = overlay[overlay["area_intersection"] > 0]

    # 5) Numerators: y_i * A_i∩n
    numerators = overlay[yield_cols].multiply(overlay["area_intersection"], axis=0)

    # 6) Sum per region
    numer_by_region = numerators.groupby(overlay["name"]).sum(min_count=1)
    denom_by_region = overlay.groupby("name")["area_intersection"].sum()

    # 7) Weighted average
    yields_regions = numer_by_region.div(denom_by_region, axis=0)

    # 8) Ensure one row per region (regions with no valid overlaps → NaN)
    yields_regions = yields_regions.reindex(regs["name"].values)

    return yields_regions


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

    # Address the missign countries (df_nuts2 contains all NUTS2 + missing shapes (with NaNs))
    missing_countries = ["AL", "RS", "BA", "XK"]
    missing_shapes = df_nuts2.loc[missing_countries, ["geometry"]]

    # Impute yields for missing shapes basest of nearest valid Nuts2
    yield_cols = [
        "Bioethanol barley, wheat, grain maize, oats, other cereals and rye",
        "Sugar from sugar beet",
        "Rape seed",
        "perennials",
    ]
    imputed_missing = impute_missing_values(df_nuts2, missing_shapes, yield_cols)
    df_nuts2 = pd.concat([
        df_nuts2.drop(index=missing_countries),
        imputed_missing
    ])

    # convert nuts2 yields to regions
    df = convert_nuts2_to_regions_yields(df_nuts2, regions)

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
