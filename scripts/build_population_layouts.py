# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build mapping between cutout grid cells and population (total, urban, rural).
"""

import logging

import atlite
import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_population_layouts",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    coco.logging.getLogger().setLevel(coco.logging.CRITICAL)

    cutout = atlite.Cutout(snakemake.input.cutout)

    grid_cells = cutout.grid.geometry

    # nuts3 has columns country, gdp, pop, geometry
    # population is given in dimensions of 1e3=k
    nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index("index")

    # Indicator matrix NUTS3 -> grid cells
    I = atlite.cutout.compute_indicatormatrix(nuts3.geometry, grid_cells)  # noqa: E741

    # Indicator matrix grid_cells -> NUTS3; inprinciple Iinv*I is identity
    # but imprecisions mean not perfect
    Iinv = cutout.indicatormatrix(nuts3.geometry)

    countries = np.sort(nuts3.country.unique())

    urban_fraction = pd.read_csv(snakemake.input.urban_percent, skiprows=4)
    iso3 = urban_fraction["Country Code"]
    urban_fraction["iso2"] = cc.convert(names=iso3, src="ISO3", to="ISO2")
    urban_fraction = (
        urban_fraction.query("iso2 in @countries").set_index("iso2")["2019"].div(100)
    )
    if "XK" in countries:
        urban_fraction["XK"] = urban_fraction["RS"]

    # population in each grid cell
    pop_cells = pd.Series(I.dot(nuts3["pop"]))

    # in km^2
    cell_areas = grid_cells.to_crs(3035).area / 1e6

    # pop per km^2
    density_cells = pop_cells / cell_areas

    # rural or urban population in grid cell
    pop_rural = pd.Series(0.0, density_cells.index)
    pop_urban = pd.Series(0.0, density_cells.index)

    for ct in countries:
        logger.debug(
            f"The urbanization rate for {ct} is {round(urban_fraction[ct] * 100)}%"
        )

        indicator_nuts3_ct = nuts3.country.apply(lambda x: 1.0 if x == ct else 0.0)

        indicator_cells_ct = pd.Series(Iinv.T.dot(indicator_nuts3_ct))

        density_cells_ct = indicator_cells_ct * density_cells

        pop_cells_ct = indicator_cells_ct * pop_cells

        # correct for imprecision of Iinv*I
        pop_ct = nuts3.loc[nuts3.country == ct, "pop"].sum()
        if pop_cells_ct.sum() != 0:
            pop_cells_ct *= pop_ct / pop_cells_ct.sum()

        # The first low density grid cells to reach rural fraction are rural
        asc_density_i = density_cells_ct.sort_values().index
        asc_density_cumsum = (
            pop_cells_ct.iloc[asc_density_i].cumsum() / pop_cells_ct.sum()
        )
        rural_fraction_ct = 1 - urban_fraction[ct]
        pop_ct_rural_b = asc_density_cumsum < rural_fraction_ct
        pop_ct_urban_b = ~pop_ct_rural_b

        pop_ct_rural_b[indicator_cells_ct == 0.0] = False
        pop_ct_urban_b[indicator_cells_ct == 0.0] = False

        pop_rural += pop_cells_ct.where(pop_ct_rural_b, 0.0)
        pop_urban += pop_cells_ct.where(pop_ct_urban_b, 0.0)

    pop_cells = {"total": pop_cells}
    pop_cells["rural"] = pop_rural
    pop_cells["urban"] = pop_urban

    for key, pop in pop_cells.items():
        ycoords = ("y", cutout.coords["y"].data)
        xcoords = ("x", cutout.coords["x"].data)
        values = pop.values.reshape(cutout.shape)
        layout = xr.DataArray(values, [ycoords, xcoords])

        layout.to_netcdf(snakemake.output[f"pop_layout_{key}"])
