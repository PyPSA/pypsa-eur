"""Cluster gas network."""

import logging
logger = logging.getLogger(__name__)

import pandas as pd
import geopandas as gpd

from shapely import wkt
from pypsa.geo import haversine_pts
from packaging.version import Version, parse


def concat_gdf(gdf_list, crs='EPSG:4326'):
    """Concatenate multiple geopandas dataframes with common coordinate reference system (crs)."""
    return gpd.GeoDataFrame(pd.concat(gdf_list), crs=crs)


def load_bus_regions(onshore_path, offshore_path):
    """Load pypsa-eur on- and offshore regions and concat."""

    bus_regions_offshore = gpd.read_file(offshore_path)
    bus_regions_onshore = gpd.read_file(onshore_path)
    bus_regions = concat_gdf([bus_regions_offshore, bus_regions_onshore])
    bus_regions = bus_regions.dissolve(by='name', aggfunc='sum')

    return bus_regions


def build_clustered_gas_network(df, bus_regions, length_factor=1.25):
    
    for i in [0,1]:

        gdf = gpd.GeoDataFrame(geometry=df[f"point{i}"], crs="EPSG:4326")

        kws = dict(op="within") if parse(gpd.__version__) < Version('0.10') else dict(predicate="within")
        bus_mapping = gpd.sjoin(gdf, bus_regions, how="left", **kws).index_right
        bus_mapping = bus_mapping.groupby(bus_mapping.index).first()

        df[f"bus{i}"] = bus_mapping

        df[f"point{i}"] = df[f"bus{i}"].map(bus_regions.to_crs(3035).centroid.to_crs(4326))

    # drop pipes where not both buses are inside regions
    df = df.loc[~df.bus0.isna() & ~df.bus1.isna()]

    # drop pipes within the same region
    df = df.loc[df.bus1 != df.bus0]

    # recalculate lengths as center to center * length factor
    df["length"] = df.apply(
        lambda p: length_factor * haversine_pts(
            [p.point0.x, p.point0.y],
            [p.point1.x, p.point1.y]
        ), axis=1
    )

    # tidy and create new numbered index
    df.drop(["point0", "point1"], axis=1, inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def reindex_pipes(df):

    def make_index(x):
        connector = " <-> " if x.bidirectional else " -> "
        return "gas pipeline " + x.bus0 + connector + x.bus1

    df.index = df.apply(make_index, axis=1)

    df["p_min_pu"] = df.bidirectional.apply(lambda bi: -1 if bi else 0)
    df.drop("bidirectional", axis=1, inplace=True)

    df.sort_index(axis=1, inplace=True)


def aggregate_parallel_pipes(df):

    strategies = {
        'bus0': 'first',
        'bus1': 'first',
        "p_nom": 'sum',
        "p_nom_diameter": 'sum',
        "max_pressure_bar": "mean",
        "build_year": "mean",
        "diameter_mm": "mean",
        "length": 'mean',
        'name': ' '.join,
        "p_min_pu": 'min',
    }
    return df.groupby(df.index).agg(strategies)


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'cluster_gas_network',
            simpl='',
            clusters='37'
        )

    logging.basicConfig(level=snakemake.config['logging_level'])

    fn = snakemake.input.cleaned_gas_network
    df = pd.read_csv(fn, index_col=0)
    for col in ["point0", "point1"]:
        df[col] = df[col].apply(wkt.loads)

    bus_regions = load_bus_regions(
        snakemake.input.regions_onshore,
        snakemake.input.regions_offshore
    )

    gas_network = build_clustered_gas_network(df, bus_regions)

    reindex_pipes(gas_network)
    gas_network = aggregate_parallel_pipes(gas_network)

    gas_network.to_csv(snakemake.output.clustered_gas_network)