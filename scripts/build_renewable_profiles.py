# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculates for each clustered region the (i) installable capacity (based on
land-use from :mod:`determine_availability_matrix`), (ii) the available
generation time series (based on weather data), and (iii) the average distance
from the node for onshore wind, AC-connected offshore wind, DC-connected
offshore wind and solar PV generators.

.. note:: Hydroelectric profiles are built in script :mod:`build_hydro_profiles`.

Outputs
-------

- ``resources/profile_{technology}.nc`` with the following structure

    ===================  ==========  =========================================================
    Field                Dimensions  Description
    ===================  ==========  =========================================================
    profile              bus, time   the per unit hourly availability factors for each bus
    -------------------  ----------  ---------------------------------------------------------
    p_nom_max            bus         maximal installable capacity at the bus (in MW)
    -------------------  ----------  ---------------------------------------------------------
    average_distance     bus         average distance of units in the region to the
                                     grid bus for onshore technologies and to the shoreline
                                     for offshore technologies (in km)
    ===================  ==========  =========================================================

    - **profile**

    .. image:: img/profile_ts.png
        :scale: 33 %
        :align: center

    - **p_nom_max**

    .. image:: img/p_nom_max_hist.png
        :scale: 33 %
        :align: center

    - **average_distance**

    .. image:: img/distance_hist.png
        :scale: 33 %
        :align: center

Description
-----------

This script functions at two main spatial resolutions: the resolution of the
clustered network regions, and the resolution of the cutout grid cells for the
weather data. Typically the weather data grid is finer than the network regions,
so we have to work out the distribution of generators across the grid cells
within each region. This is done by taking account of a combination of the
available land at each grid cell (computed in
:mod:`determine_availability_matrix`) and the capacity factor there.

Based on the availability matrix, the script first computes how much of the
technology can be installed at each cutout grid cell. To compute the layout of
generators in each clustered region, the installable potential in each grid cell
is multiplied with the capacity factor at each grid cell. This is done since we
assume more generators are installed at cells with a higher capacity factor.

.. image:: img/offwinddc-gridcell.png
    :scale: 50 %
    :align: center

.. image:: img/offwindac-gridcell.png
    :scale: 50 %
    :align: center

.. image:: img/onwind-gridcell.png
    :scale: 50 %
    :align: center

.. image:: img/solar-gridcell.png
    :scale: 50 %
    :align: center

This layout is then used to compute the generation availability time series from
the weather data cutout from ``atlite``.

The maximal installable potential for the node (`p_nom_max`) is computed by
adding up the installable potentials of the individual grid cells.
"""

import logging
import time

import atlite
import geopandas as gpd
import xarray as xr
from _helpers import configure_logging, get_snapshots, set_scenario_config
from build_shapes import _simplify_polys
from dask.distributed import Client

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_renewable_profiles", clusters=38, technology="offwind-ac"
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    noprogress = snakemake.config["run"].get("disable_progressbar", True)
    noprogress = noprogress or not snakemake.config["atlite"]["show_progress"]
    technology = snakemake.wildcards.technology
    params = snakemake.params.renewable[technology]
    resource = params["resource"]  # pv panel params / wind turbine params
    resource["show_progress"] = not noprogress

    tech = next(t for t in ["panel", "turbine"] if t in resource)
    models = resource[tech]
    if not isinstance(models, dict):
        models = {0: models}
    resource[tech] = models[next(iter(models))]

    correction_factor = params.get("correction_factor", 1.0)
    capacity_per_sqkm = params["capacity_per_sqkm"]

    if correction_factor != 1.0:
        logger.info(f"correction_factor is set as {correction_factor}")

    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
    else:
        client = None

    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)

    cutout = atlite.Cutout(snakemake.input.cutout).sel(time=sns)

    availability = xr.open_dataarray(snakemake.input.availability_matrix)

    regions = gpd.read_file(snakemake.input.regions)
    assert not regions.empty, (
        f"List of regions in {snakemake.input.regions} is empty, please "
        "disable the corresponding renewable technology"
    )
    # do not pull up, set_index does not work if geo dataframe is empty
    regions = regions.set_index("name").rename_axis("bus")
    if snakemake.wildcards.technology.startswith("offwind"):
        # for offshore regions, the shortest distance to the shoreline is used
        offshore_regions = availability.coords["bus"].values
        regions = regions.loc[offshore_regions]
        regions = regions.map(lambda g: _simplify_polys(g, minarea=1)).set_crs(
            regions.crs
        )
    else:
        # for onshore regions, the representative point of the region is used
        regions = regions.representative_point()
    regions = regions.geometry.to_crs(3035)
    buses = regions.index

    area = cutout.grid.to_crs(3035).area / 1e6
    area = xr.DataArray(
        area.values.reshape(cutout.shape), [cutout.coords["y"], cutout.coords["x"]]
    )

    func = getattr(cutout, resource.pop("method"))
    if client is not None:
        resource["dask_kwargs"] = {"scheduler": client}

    logger.info(f"Calculate average capacity factor for technology {technology}...")
    start = time.time()

    capacity_factor = correction_factor * func(capacity_factor=True, **resource)
    layout = capacity_factor * area * capacity_per_sqkm

    duration = time.time() - start
    logger.info(
        f"Completed average capacity factor calculation for technology {technology} ({duration:2.2f}s)"
    )

    profiles = []
    for year, model in models.items():
        logger.info(
            f"Calculate weighted capacity factor time series for model {model} for technology {technology}..."
        )
        start = time.time()

        resource[tech] = model

        profile = func(
            matrix=availability.stack(spatial=["y", "x"]),
            layout=layout,
            index=buses,
            per_unit=True,
            return_capacity=False,
            **resource,
        )

        dim = {"year": [year]}
        profile = profile.expand_dims(dim)

        profiles.append(profile.rename("profile"))

        duration = time.time() - start
        logger.info(
            f"Completed weighted capacity factor time series calculation for model {model} for technology {technology} ({duration:2.2f}s)"
        )

    profiles = xr.merge(profiles)

    logger.info(f"Calculating maximal capacity per bus for technology {technology}")
    p_nom_max = capacity_per_sqkm * availability @ area

    logger.info(f"Calculate average distances for technology {technology}.")
    layoutmatrix = (layout * availability).stack(spatial=["y", "x"])

    coords = cutout.grid.representative_point().to_crs(3035)

    average_distance = []
    for bus in buses:
        row = layoutmatrix.sel(bus=bus).data
        nz_b = row != 0
        row = row[nz_b]
        co = coords[nz_b]
        distances = co.distance(regions[bus]).div(1e3)  # km
        average_distance.append((distances * (row / row.sum())).sum())

    average_distance = xr.DataArray(average_distance, [buses])

    ds = xr.merge(
        [
            correction_factor * profiles,
            p_nom_max.rename("p_nom_max"),
            average_distance.rename("average_distance"),
        ]
    )
    # select only buses with some capacity and minimal capacity factor
    mean_profile = ds["profile"].mean("time")
    if "year" in ds.indexes:
        mean_profile = mean_profile.max("year")

    ds = ds.sel(
        bus=(
            (mean_profile > params.get("min_p_max_pu", 0.0))
            & (ds["p_nom_max"] > params.get("min_p_nom_max", 0.0))
        )
    )

    if "clip_p_max_pu" in params:
        min_p_max_pu = params["clip_p_max_pu"]
        ds["profile"] = ds["profile"].where(ds["profile"] >= min_p_max_pu, 0)

    ds.to_netcdf(snakemake.output.profile)

    if client is not None:
        client.shutdown()
