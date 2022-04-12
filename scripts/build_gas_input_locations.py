"""
Build import locations for fossil gas from entry-points, LNG terminals and production sites.
"""

import logging
logger = logging.getLogger(__name__)

import pandas as pd
import geopandas as gpd
from shapely import wkt

from cluster_gas_network import load_bus_regions


def read_scigrid_gas(fn):
    df = gpd.read_file(fn)
    df = pd.concat([df, df.param.apply(pd.Series)], axis=1)
    df.drop(["param", "uncertainty", "method"], axis=1, inplace=True)
    return df


def build_gas_input_locations(lng_fn, planned_lng_fn, entry_fn, prod_fn, countries):
    
    # LNG terminals
    lng = read_scigrid_gas(lng_fn)
    planned_lng = pd.read_csv(planned_lng_fn)
    planned_lng.geometry = planned_lng.geometry.apply(wkt.loads)
    planned_lng = gpd.GeoDataFrame(planned_lng, crs=4326)
    lng = pd.concat([lng, planned_lng], ignore_index=True)

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

    conversion_factor = 437.5 # MCM/day to MWh/h
    lng["p_nom"] = lng["max_cap_store2pipe_M_m3_per_d"] * conversion_factor
    entry["p_nom"] = entry["max_cap_from_to_M_m3_per_d"] * conversion_factor
    prod["p_nom"] = prod["max_supply_M_m3_per_d"] * conversion_factor

    lng["type"] = "lng"
    entry["type"] = "pipeline"
    prod["type"] = "production"

    sel = ["geometry", "p_nom", "type"]

    return pd.concat([prod[sel], entry[sel], lng[sel]], ignore_index=True)


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_gas_import_locations',
            simpl='',
            clusters='37',
        )

    logging.basicConfig(level=snakemake.config['logging_level'])

    regions = load_bus_regions(
        snakemake.input.regions_onshore,
        snakemake.input.regions_offshore
    )

    # add a buffer to eastern countries because some
    # entry points are still in Russian or Ukrainian territory.
    buffer = 9000 # meters
    eastern_countries = ['FI', 'EE', 'LT', 'LV', 'PL', 'SK', 'HU', 'RO']
    add_buffer_b = regions.index.str[:2].isin(eastern_countries)
    regions.loc[add_buffer_b] = regions[add_buffer_b].to_crs(3035).buffer(buffer).to_crs(4326)

    countries = regions.index.str[:2].unique().str.replace("GB", "UK")

    gas_input_locations = build_gas_input_locations(
        snakemake.input.lng,
        snakemake.input.planned_lng,
        snakemake.input.entry,
        snakemake.input.production,
        countries
    )

    gas_input_nodes = gpd.sjoin(gas_input_locations, regions, how='left')

    gas_input_nodes.rename(columns={"index_right": "bus"}, inplace=True)

    gas_input_nodes.to_file(snakemake.output.gas_input_nodes, driver='GeoJSON')

    gas_input_nodes_s = gas_input_nodes.groupby(["bus", "type"])["p_nom"].sum().unstack()
    gas_input_nodes_s.columns.name = "p_nom"

    gas_input_nodes_s.to_csv(snakemake.output.gas_input_nodes_simplified)