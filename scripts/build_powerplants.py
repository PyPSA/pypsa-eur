# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


"""
Retrieves conventional powerplant capacities and locations from
`powerplantmatching <https://github.com/PyPSA/powerplantmatching>`_, assigns
these to buses and creates a ``.csv`` file. It is possible to amend the
powerplant database with custom entries provided in
``data/custom_powerplants.csv``.
Lastly, for every substation, powerplants with zero-initial capacity can be added for certain fuel types automatically.

Outputs
-------

- ``resource/powerplants_s_{clusters}.csv``: A list of conventional power plants (i.e. neither wind nor solar) with fields for name, fuel type, technology, country, capacity in MW, duration, commissioning year, retrofit year, latitude, longitude, and dam information as documented in the `powerplantmatching README <https://github.com/PyPSA/powerplantmatching/blob/master/README.md>`_; additionally it includes information on the closest substation/bus in ``networks/base_s_{clusters}.nc``.

    .. image:: img/powerplantmatching.png
        :scale: 30 %

    **Source:** `powerplantmatching on GitHub <https://github.com/PyPSA/powerplantmatching>`_

Description
-----------

The configuration options ``electricity: powerplants_filter`` and ``electricity: custom_powerplants`` can be used to control whether data should be retrieved from the original powerplants database or from custom amendments. These specify `pandas.query <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.query.html>`_ commands.
In addition the configuration option ``electricity: everywhere_powerplants`` can be used to place powerplants with zero-initial capacity of certain fuel types at all substations.

1. Adding all powerplants from custom:

    .. code:: yaml

        powerplants_filter: false
        custom_powerplants: true

2. Replacing powerplants in e.g. Germany by custom data:

    .. code:: yaml

        powerplants_filter: Country not in ['Germany']
        custom_powerplants: true

    or

    .. code:: yaml

        powerplants_filter: Country not in ['Germany']
        custom_powerplants: Country in ['Germany']


3. Adding additional built year constraints:

    .. code:: yaml

        powerplants_filter: Country not in ['Germany'] and YearCommissioned <= 2015
        custom_powerplants: YearCommissioned <= 2015

4. Adding powerplants at all substations for 4 conventional carrier types:

    .. code:: yaml

        everywhere_powerplants: ['Natural Gas', 'Coal', 'nuclear', 'OCGT']
"""

import itertools
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
from shapely.geometry import MultiPolygon, Polygon

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def add_custom_powerplants(ppl, custom_powerplants, custom_ppl_query=False):
    if not custom_ppl_query:
        return ppl
    add_ppls = pd.read_csv(custom_powerplants, dtype={"bus": "str"})
    if isinstance(custom_ppl_query, str):
        add_ppls.query(custom_ppl_query, inplace=True)
    return pd.concat(
        [ppl, add_ppls], sort=False, ignore_index=True, verify_integrity=True
    )


def add_everywhere_powerplants(ppl, substations, everywhere_powerplants):
    # Create a dataframe with "everywhere_powerplants" of stated carriers at the location of all substations
    everywhere_ppl = (
        pd.DataFrame(
            itertools.product(substations.index.values, everywhere_powerplants),
            columns=["substation_index", "Fueltype"],
        ).merge(
            substations[["x", "y", "country"]],
            left_on="substation_index",
            right_index=True,
        )
    ).drop(columns="substation_index")

    # PPL uses different columns names compared to substations dataframe -> rename
    everywhere_ppl = everywhere_ppl.rename(
        columns={"x": "lon", "y": "lat", "country": "Country"}
    )

    # Add default values for the powerplants
    everywhere_ppl["Name"] = (
        "Automatically added everywhere-powerplant " + everywhere_ppl.Fueltype
    )
    everywhere_ppl["Set"] = "PP"
    everywhere_ppl["Technology"] = everywhere_ppl["Fueltype"]
    everywhere_ppl["Capacity"] = 0.0

    # Assign plausible values for the commissioning and decommissioning years
    # required for multi-year models
    everywhere_ppl["DateIn"] = ppl["DateIn"].min()
    everywhere_ppl["DateOut"] = ppl["DateOut"].max()

    # NaN values for efficiency will be replaced by the generic efficiency by attach_conventional_generators(...) in add_electricity.py later
    everywhere_ppl["Efficiency"] = np.nan

    return pd.concat(
        [ppl, everywhere_ppl], sort=False, ignore_index=True, verify_integrity=True
    )


def replace_natural_gas_technology(df):
    mapping = {
        "Steam Turbine": "CCGT",
        "Combustion Engine": "OCGT",
        "Not Found": "CCGT",
    }
    tech = df.Technology.replace(mapping).fillna("CCGT")
    return df.Technology.mask(df.Fueltype == "Natural Gas", tech)


def replace_natural_gas_fueltype(df: pd.DataFrame) -> pd.Series:
    return df.Fueltype.mask(
        (df.Technology == "OCGT") | (df.Technology == "CCGT"), "Natural Gas"
    )


def fill_unoccupied_holes(gdf: gpd.GeoDataFrame) -> gpd.GeoSeries:
    def _fill_poly(poly, idx):
        if not poly.interiors:
            return poly
        kept = [h for h in poly.interiors if gdf.drop(idx).intersects(Polygon(h)).any()]
        return Polygon(poly.exterior, kept)

    result = gdf.geometry.copy()
    for idx in gdf.index:
        g = gdf.geometry[idx]
        if g.geom_type == "Polygon":
            result[idx] = _fill_poly(g, idx)
        elif g.geom_type == "MultiPolygon":
            result[idx] = MultiPolygon([_fill_poly(p, idx) for p in g.geoms])
    return result


def map_to_country_bus(
    ppl: gpd.GeoDataFrame, regions: gpd.GeoDataFrame, max_distance: float = 10000
) -> gpd.GeoDataFrame:
    """
    Assign power plants to region buses of the same country.

    First, spatial join is performed per country to avoid cross-border
    misassignment. Remaining unmatched plants are assigned via nearest
    neighbor (max 10000m) within the same country.
    """
    assigned = []
    unmatched = []

    for country, plants in ppl.groupby("Country"):
        country_regions = regions[regions.index.str[:2] == country]
        joined = (
            plants.sjoin(country_regions)
            .rename(columns={"name": "bus"})
            .reindex(plants.index)
        )
        assigned.append(joined.dropna(subset=["bus"]))
        missing = joined[joined["bus"].isna()]
        if not missing.empty:
            unmatched.append(plants.loc[missing.index])

    if unmatched:
        unmatched = pd.concat(unmatched)
        for country, plants in unmatched.groupby("Country"):
            country_regions = regions[regions.index.str[:2] == country]
            nearest = (
                plants.to_crs(3035)
                .sjoin_nearest(country_regions.to_crs(3035), max_distance=max_distance)
                .rename(columns={"name": "bus"})
                .to_crs(4326)
            )
            missing = plants.index.difference(nearest.index)
            print(country, missing)
            nearest = pd.concat([nearest, plants.loc[missing]])
            assigned.append(nearest)

    return pd.concat(assigned)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_powerplants", clusters=256)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    countries = snakemake.params.countries

    fn_onshore = snakemake.input.regions_onshore
    fn_offshore = snakemake.input.regions_offshore

    regions = pd.concat([gpd.read_file(fn_onshore), gpd.read_file(fn_offshore)])
    regions = regions.dissolve("name")
    regions["geometry"] = fill_unoccupied_holes(regions)

    # Steps copied from PPM: Usually run by PPM when using pm.powerplants(...) from cache
    ppl = (
        pd.read_csv(snakemake.input.powerplants, index_col=0, header=[0])
        .pipe(pm.collection.parse_string_to_dict, ["projectID", "EIC"])
        .pipe(pm.collection.set_column_name, "Matched Data")
    )
    ppl = (
        ppl.powerplant.convert_country_to_alpha2()
        .query("Country in @countries")
        .assign(Technology=replace_natural_gas_technology)
        .assign(Fueltype=replace_natural_gas_fueltype)
        .replace({"Solid Biomass": "Bioenergy", "Biogas": "Bioenergy"})
    )

    ppl_query = snakemake.params.powerplants_filter
    if isinstance(ppl_query, str):
        ppl.query(ppl_query, inplace=True)

    # add carriers from own powerplant files:
    custom_ppl_query = snakemake.params.custom_powerplants
    ppl = add_custom_powerplants(
        ppl, snakemake.input.custom_powerplants, custom_ppl_query
    )

    if countries_wo_ppl := set(countries) - set(ppl.Country.unique()):
        logger.warning(f"No powerplants known in: {', '.join(countries_wo_ppl)}")

    # Add "everywhere powerplants" to all bus locations
    ppl = add_everywhere_powerplants(
        ppl, n.buses, snakemake.params.everywhere_powerplants
    )

    ppl = ppl.dropna(subset=["lat", "lon"])

    ppl = gpd.GeoDataFrame(ppl, geometry=gpd.points_from_xy(ppl.lon, ppl.lat), crs=4326)

    ppl = map_to_country_bus(ppl, regions)

    bus_null_b = ppl["bus"].isnull()
    if bus_null_b.any():
        stats = (
            ppl.loc[bus_null_b]
            .groupby(by=["Country", "Fueltype"])
            .Capacity.sum()
            .sort_values(ascending=False)
        )
        logger.warning(
            f"Couldn't assign sufficiently close region for {bus_null_b.sum()} powerplants.\n"
            f"Removing the following capacities (MW) from the powerplants dataset:\n {stats}"
        )
        ppl = ppl[~bus_null_b]

    ppl.reset_index(drop=True).to_csv(snakemake.output[0])
