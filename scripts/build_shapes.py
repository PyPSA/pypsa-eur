# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates GIS shape files of the countries, exclusive economic zones and `NUTS3 <
https://en.wikipedia.org/wiki/Nomenclature_of_Territorial_Units_for_Statistics>
`_ areas.

Relevant Settings
-----------------

.. code:: yaml

    countries:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`toplevel_cf`

Inputs
------

- ``data/bundle/naturalearth/ne_10m_admin_0_countries.shp``: World country shapes

    .. image:: img/countries.png
        :scale: 33 %

- ``data/eez/World_EEZ_v12_20231025_gpkg/eez_v12.gpkg   ``: World `exclusive economic zones <https://en.wikipedia.org/wiki/Exclusive_economic_zone>`_ (EEZ)

    .. image:: img/eez.png
        :scale: 33 %

- ``data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp``: Europe NUTS3 regions

    .. image:: img/nuts3.png
        :scale: 33 %

- ``data/bundle/nama_10r_3popgdp.tsv.gz``: Average annual population by NUTS3 region (`eurostat <http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nama_10r_3popgdp&lang=en>`__)
- ``data/bundle/nama_10r_3gdp.tsv.gz``: Gross domestic product (GDP) by NUTS 3 regions (`eurostat <http://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=nama_10r_3gdp&lang=en>`__)
- ``data/ch_cantons.csv``: Mapping between Swiss Cantons and NUTS3 regions
- ``data/bundle/je-e-21.03.02.xls``: Population and GDP data per Canton (`BFS - Swiss Federal Statistical Office <https://www.bfs.admin.ch/bfs/en/home/news/whats-new.assetdetail.7786557.html>`_ )

Outputs
-------

- ``resources/country_shapes.geojson``: country shapes out of country selection

    .. image:: img/country_shapes.png
        :scale: 33 %

- ``resources/offshore_shapes.geojson``: EEZ shapes out of country selection

    .. image:: img/offshore_shapes.png
        :scale: 33 %

- ``resources/europe_shape.geojson``: Shape of Europe including countries and EEZ

    .. image:: img/europe_shape.png
        :scale: 33 %

- ``resources/nuts3_shapes.geojson``: NUTS3 shapes out of country selection including population and GDP data.

    .. image:: img/nuts3_shapes.png
        :scale: 33 %

Description
-----------
"""
from itertools import takewhile
from operator import attrgetter
from rasterio.mask import mask
from shapely.geometry import box, MultiPolygon, Polygon

import country_converter as coco
import geopandas as gpd
import logging
import numpy as np
import pandas as pd
import rasterio
import unicodedata
import xarray as xr

from _helpers import configure_logging, set_scenario_config


logger = logging.getLogger(__name__)
cc = coco.CountryConverter()


GDP_YEAR=2019
POP_YEAR=2019
DROP_REGIONS = [
    "ES703",
    "ES704",
    "ES705",
    "ES706",
    "ES707",
    "ES708",
    "ES709",
    "ES630",
    "ES640",
    "FRY10",
    "FRY20",
    "FRY30",
    "FRY40",
    "FRY50",
    "NO0B1",
    "NO0B2",
    "PT200",
    "PT300",
]
OTHER_GDP_TOTAL_2019 = {    # in bn. USD
    "BA": 20.48,            # World Bank
    "MD": 11.74,            # World Bank
    "UA": 153.9,            # World Bank
    "XK": 7.9,              # https://de.statista.com/statistik/daten/studie/415738/umfrage/bruttoinlandsprodukt-bip-des-kosovo/
}
OTHER_POP_2019 = {          # in 1000 persons
    "BA": 3361,             # World Bank
    "MD": 2664,             # World Bank
    "UA": 44470,            # World Bank
    "XK": 1782,             # World Bank
}
EXCHANGE_EUR_USD_2019 = 1.1


def _simplify_polys(polys, minarea=0.1, tolerance=None, filterremote=True):
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys.geoms, key=attrgetter("area"), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area / (2.0 * np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon(
                [
                    p
                    for p in takewhile(lambda p: p.area > minarea, polys)
                    if not filterremote or (mainpoly.distance(p) < mainlength)
                ]
            )
        else:
            polys = mainpoly
    if tolerance is not None:
        polys = polys.simplify(tolerance=tolerance)
    return polys


def eez(eez, country_list):
    df = gpd.read_file(eez)
    iso3_list = cc.convert(country_list, src="ISO2", to="ISO3")
    pol_type = ["200NM", "Overlapping claim"]
    df = df.query("ISO_TER1 in @iso3_list and POL_TYPE in @pol_type").copy()
    df["name"] = cc.convert(df.ISO_TER1, src="ISO3", to="ISO2")
    s = df.set_index("name").geometry.map(
        lambda s: _simplify_polys(s, filterremote=False)
    )
    s = s.to_frame("geometry").set_crs(df.crs)
    s.index.name = "name"
    return s


def country_cover(country_shapes, eez_shapes=None):
    shapes = country_shapes
    if eez_shapes is not None:
        shapes = pd.concat([shapes, eez_shapes])
    europe_shape = shapes.union_all()
    if isinstance(europe_shape, MultiPolygon):
        europe_shape = max(europe_shape.geoms, key=attrgetter("area"))
    return Polygon(shell=europe_shape.exterior)


def normalise_text(text):
    """
    Removes diacritics from non-standard Latin letters, converts them to their
    closest standard 26-letter Latin equivalents, and removes asterisks (*) from the text.

    Args:
        text (str): Input string to normalize.

    Returns:
        str: Normalized string with only standard Latin letters and no asterisks.
    """
    # Normalize Unicode to decompose characters (e.g., č -> c + ̌)
    text = unicodedata.normalize('NFD', text)
    # Remove diacritical marks by filtering out characters of the 'Mn' (Mark, Nonspacing) category
    text = ''.join(char for char in text if unicodedata.category(char) != 'Mn')
    # Remove asterisks
    text = text.replace('*', '')
    # Optionally, ensure only ASCII characters remain
    text = ''.join(char for char in text if char.isascii())
    return text


def create_regions(
    country_list,
    nuts3_path,
    ba_adm1_path,
    md_adm1_path,
    ua_adm1_path,
    xk_adm1_path,
    offshore_shapes,
):
    """
    Create regions by processing NUTS and non-NUTS geographical shapes.

    Parameters:
        - country_list (list): List of country codes to include.
        - nuts3_path (str): Path to the NUTS3 2021 shapefile.
        - ba_adm1_path (str): Path to adm1 boundaries for Bosnia and  Herzegovina.
        - md_adm1_path (str): Path to adm1 boundaries for Moldova.
        - ua_adm1_path (str): Path to adm1 boundaries for Ukraine.
        - xk_adm1_path (str): Path to adm1 boundaries for Kosovo.
        - offshore_shapes (geopandas.GeoDataFrame): Geographical shapes of the exclusive economic zones.
    
    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame containing the processed regions with columns:
            - id: Region identifier.
            - country: Country code.
            - name: Region name.
            - geometry: Geometrical shape of the region.
            - level1: Level 1 region identifier.
            - level2: Level 2 region identifier.
            - level3: Level 3 region identifier.
    """
    # Prepare NUTS shapes
    logger.info("Processing NUTS regions.")
    regions = gpd.read_file(nuts3_path)
    regions.loc[regions.CNTR_CODE == "EL", "CNTR_CODE"] = "GR" # Rename "EL" to "GR
    regions["NUTS_ID"] = regions["NUTS_ID"].str.replace("EL", "GR")
    regions.loc[regions.CNTR_CODE == "UK", "CNTR_CODE"] = "GB" # Rename "UK" to "GB"
    regions["NUTS_ID"] = regions["NUTS_ID"].str.replace("UK", "GB")

    # Only include countries in the config
    regions = regions.query("CNTR_CODE in @country_list")

    # Create new df
    regions = regions[["NUTS_ID", "CNTR_CODE", "NAME_LATN", "geometry"]]
     
    # Rename columns and add level columns
    regions = regions.rename(columns={"NUTS_ID": "id", "CNTR_CODE": "country", "NAME_LATN": "name"})
    
    # Normalise text
    regions["id"] = regions["id"].apply(normalise_text)
    
    regions["level1"] = regions["id"].str[:3]
    regions["level2"] = regions["id"].str[:4]
    regions["level3"] = regions["id"]

    # Non NUTS countries
    logger.info("Processing non-NUTS regions.")
    
    ba_adm1 = gpd.read_file(ba_adm1_path)
    md_adm1 = gpd.read_file(md_adm1_path)
    ua_adm1 = gpd.read_file(ua_adm1_path)
    xk_adm1 = gpd.read_file(xk_adm1_path)

    regions_non_nuts = pd.concat([ba_adm1, md_adm1, ua_adm1, xk_adm1])
    regions_non_nuts = regions_non_nuts.drop(columns=["osm_id"])

    # Normalise text
    regions_non_nuts["id"] = regions_non_nuts["id"].apply(normalise_text)
    regions_non_nuts["name"] = regions_non_nuts["name"].apply(normalise_text)

    # Add level columns
    regions_non_nuts["level1"] = regions_non_nuts["id"]
    regions_non_nuts["level2"] = regions_non_nuts["id"]
    regions_non_nuts["level3"] = regions_non_nuts["id"]

    # Clip regions by non-NUTS shapes
    regions["geometry"] = regions["geometry"].difference(regions_non_nuts.geometry.union_all())

    # Concatenate NUTS and non-NUTS regions
    logger.info("Harmonising NUTS and non-NUTS regions.")
    regions = pd.concat([regions, regions_non_nuts])
    regions.set_index("id", inplace=True)

    # Drop regions out of geographical scope
    regions = regions.drop(DROP_REGIONS, errors="ignore")

    # Clip regions by offshore shapes
    logger.info("Clipping regions by offshore shapes.")
    regions["geometry"] = regions["geometry"].difference(offshore_shapes.geometry.union_all())

    return regions


def calc_gdp_pop(country, regions, gdp_non_nuts3, pop_non_nuts3):
    """
    Calculate the GDP p.c. and population values for non NUTS3 regions.

    Parameters:
        - country (str): The two-letter country code of the non-NUTS3 region.
        - regions (GeoDataFrame): A GeoDataFrame containing the regions.
        - gdp_non_nuts3 (str): The file path to the dataset containing the GDP p.c values
          for non NUTS3 countries (e.g. MD, UA)
        - pop_non_nuts3 (str): The file path to the dataset containing the POP values
        for non NUTS3 countries (e.g. MD, UA)

    Returns:
    tuple: A tuple containing two GeoDataFrames:
        - gdp: A GeoDataFrame with the mean GDP p.c. values mapped to each bus.
        - pop: A GeoDataFrame with the summed POP values mapped to each bus.
    """
    region = regions.loc[regions.country == country, ["geometry"]]
    # Create a bounding box for UA, MD from region shape, including a buffer of 10000 metres
    bounding_box = (
        gpd.GeoDataFrame(geometry=[box(*region.total_bounds)], crs=region.crs)
        .to_crs(epsg=3857)
        .buffer(10000)
        .to_crs(region.crs)
    )

    # GDP Mapping
    logger.info(f"Mapping mean GDP p.c. to non-NUTS3 region: {country}")
    with xr.open_dataset(gdp_non_nuts3) as src_gdp:
        src_gdp = src_gdp.where(
            (src_gdp.longitude >= bounding_box.bounds.minx.min())
            & (src_gdp.longitude <= bounding_box.bounds.maxx.max())
            & (src_gdp.latitude >= bounding_box.bounds.miny.min())
            & (src_gdp.latitude <= bounding_box.bounds.maxy.max()),
            drop=True,
        )
        gdp = src_gdp.to_dataframe().reset_index()
    gdp = gdp.rename(columns={"GDP_per_capita_PPP": "gdp"})
    gdp = gdp[gdp.time == gdp.time.max()]
    gdp_raster = gpd.GeoDataFrame(
        gdp,
        geometry=gpd.points_from_xy(gdp.longitude, gdp.latitude),
        crs="EPSG:4326",
    )
    gdp_mapped = gpd.sjoin(gdp_raster, region, predicate="within")
    gdp = (
        gdp_mapped
        .groupby(["id"])
        .agg({"gdp": "mean"})
    )

    # Population Mapping
    logger.info(f"Mapping summed population to non-NUTS3 region: {country}")
    with rasterio.open(pop_non_nuts3) as src_pop:
        # Mask the raster with the bounding box
        out_image, out_transform = mask(src_pop, bounding_box, crop=True)
        out_meta = src_pop.meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            }
        )
    masked_data = out_image[0]  # Use the first band (rest is empty)
    row_indices, col_indices = np.where(masked_data != src_pop.nodata)
    values = masked_data[row_indices, col_indices]

    # Affine transformation from pixel coordinates to geo coordinates
    x_coords, y_coords = rasterio.transform.xy(out_transform, row_indices, col_indices)
    pop_raster = pd.DataFrame({"x": x_coords, "y": y_coords, "pop": values})
    pop_raster = gpd.GeoDataFrame(
        pop_raster,
        geometry=gpd.points_from_xy(pop_raster.x, pop_raster.y),
        crs=src_pop.crs,
    )
    pop_mapped = gpd.sjoin(pop_raster, region, predicate="within")
    pop = (
        pop_mapped.groupby(["id"])
        .agg({"pop": "sum"})
        .reset_index()
        .set_index("id")
    )
    gdp_pop = region.join(gdp).join(pop).drop(columns="geometry")
    gdp_pop.fillna(0, inplace=True)

    # Clean and rescale data to 2019 historical values
    gdp_pop["gdp"] = gdp_pop["gdp"].round(0)
    gdp_pop["pop"] = gdp_pop["pop"].div(1e3).round(0)

    gdp_pop["pop"] = (
        gdp_pop["pop"].div(gdp_pop["pop"].sum()).mul(OTHER_POP_2019[country])
        .round(0)
    )

    gdp_pop["gdp"] = (
        gdp_pop["gdp"].mul(1e9).div(gdp_pop["gdp"].sum())
        .mul(OTHER_GDP_TOTAL_2019[country])
        .div(EXCHANGE_EUR_USD_2019) /
        (1e3*gdp_pop["pop"])
    ).round(0)
        
    return gdp_pop


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_shapes")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Offshore regions
    offshore_shapes = eez(snakemake.input.eez, snakemake.params.countries)
    offshore_shapes.reset_index().to_file(snakemake.output.offshore_shapes)

    # Onshore regions
    regions = create_regions(
        snakemake.params.countries,
        snakemake.input.nuts3_2021,
        snakemake.input.ba_adm1,
        snakemake.input.md_adm1,
        snakemake.input.ua_adm1,
        snakemake.input.xk_adm1,
        offshore_shapes,
    )

    country_shapes = regions.groupby("country")["geometry"].apply(lambda x: x.union_all())
    country_shapes.crs = regions.crs
    country_shapes.index.name = "name"
    country_shapes.reset_index().to_file(snakemake.output.country_shapes)

    europe_shape = gpd.GeoDataFrame(
        geometry=[country_cover(country_shapes, offshore_shapes.geometry)],
        crs=country_shapes.crs,
    )
    europe_shape.reset_index().to_file(snakemake.output.europe_shape)

    # GDP and POP for NUTS3 regions
    # GDP
    logger.info(f"Importing JRC ARDECO GDP data for year {GDP_YEAR}.")
    nuts3_gdp = pd.read_csv(snakemake.input.nuts3_gdp, index_col=[0])
    nuts3_gdp = nuts3_gdp.query("LEVEL_ID == 3 and UNIT == 'EUR'")
    nuts3_gdp.index = nuts3_gdp.index.str.replace("UK", "GB").str.replace("EL", "GR")
    nuts3_gdp = nuts3_gdp[str(GDP_YEAR)]
    regions["gdp"] = nuts3_gdp

    # Population
    logger.info(f"Importing JRC ARDECO population data for year {POP_YEAR}.")
    nuts3_pop = pd.read_csv(snakemake.input.nuts3_pop, index_col=[0])
    nuts3_pop = nuts3_pop.query("LEVEL_ID == 3")
    nuts3_pop.index = nuts3_pop.index.str.replace("UK", "GB").str.replace("EL", "GR")
    nuts3_pop = nuts3_pop[str(POP_YEAR)]
    regions["pop"] = nuts3_pop.div(1e3).round(0)

    # GDP and POP for non-NUTS3 regions
    other_countries = {"BA", "MD", "UA", "XK"}
    other_gdp = snakemake.input.other_gdp
    other_pop = snakemake.input.other_pop

    gdp_pop = pd.concat(
        [
            calc_gdp_pop(country, regions, other_gdp, other_pop)
            for country in other_countries
        ],
        axis=0,
    )

    # Merge NUTS3 and non-NUTS3 regions
    regions.loc[gdp_pop.index, ["gdp", "pop"]] = gdp_pop[["gdp", "pop"]]

    # Resort columns and rename index
    regions = regions[["name", "level1", "level2", "level3", "gdp", "pop", "country", "geometry"]]
    regions.index.name = "index"

    # Export regions including GDP and POP data
    logger.info(f"Exporting NUTS3 and ADM1 shapes with GDP and POP values to {snakemake.output.nuts3_shapes}.")
    regions.reset_index().to_file(snakemake.output.nuts3_shapes)