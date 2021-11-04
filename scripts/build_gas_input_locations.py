"""
Build import locations for fossil gas from entry-points, LNG terminals and production sites.
"""

import logging
logger = logging.getLogger(__name__)

import pandas as pd
import geopandas as gpd


def read_scigrid_gas(fn):
    df = gpd.read_file(fn)
    df = pd.concat([df, df.param.apply(pd.Series)], axis=1)
    df.drop(["param", "uncertainty", "method"], axis=1, inplace=True)
    return df


def build_gas_input_locations(lng_fn, entry_fn, prod_fn, countries):
    
    # LNG terminals
    lng = read_scigrid_gas(lng_fn)

    # Entry points from outside the model scope
    entry = read_scigrid_gas(entry_fn)
    entry["from_country"] = entry.from_country.str.rstrip()
    entry = entry.loc[
        ~(entry.from_country.isin(countries) & entry.to_country.isin(countries)) &  # only take non-EU entries
        ~entry.name.str.contains("Tegelen") |  # malformed datapoint
        (entry.from_country == "NO")  # entries from NO to GB
    ]
    
    # production sites inside the model scope
    prod = read_scigrid_gas(prod_fn)
    prod = prod.loc[
        (prod.geometry.y > 35) &
        (prod.geometry.x < 30) &
        (prod.country_code != "DE")
    ]

    return gpd.GeoDataFrame(
        geometry=pd.concat([prod.geometry, entry.geometry, lng.geometry]).reset_index(drop=True),
        crs=4326
    )


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_gas_import_locations',
            simpl='',
            clusters='37',
        )

    logging.basicConfig(level=snakemake.config['logging_level'])

    onshore_regions = gpd.read_file(snakemake.input.regions_onshore).set_index('name')

    countries = onshore_regions.index.str[:2].unique().str.replace("GB", "UK")

    gas_input_locations = build_gas_input_locations(
        snakemake.input.lng,
        snakemake.input.entry,
        snakemake.input.production,
        countries
    )

    # recommended to use projected CRS rather than geographic CRS
    gas_input_nodes = gpd.sjoin_nearest(
        gas_input_locations.to_crs(3035),
        onshore_regions.to_crs(3035),
        how='left'
    ).index_right.unique()

    pd.Series(gas_input_nodes, name='gas_input_nodes').to_csv(snakemake.output.gas_input_nodes)