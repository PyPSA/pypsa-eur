#!/usr/bin/env python

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Calculates for each network node the
(i) installable capacity (based on land-use), (ii) the available generation time
series (based on weather data), and (iii) the average distance from the node for
onshore wind, AC-connected offshore wind, DC-connected offshore wind and solar
PV generators. In addition for offshore wind it calculates the fraction of the
grid connection which is under water.

.. note:: Hydroelectric profiles are built in script :mod:`build_hydro_profiles`.

Relevant settings
-----------------

.. code:: yaml

    snapshots:

    atlite:
        nprocesses:

    renewable:
        {technology}:
            cutout:
            corine:
            grid_codes:
            distance:
            natura:
            max_depth:
            max_shore_distance:
            min_shore_distance:
            capacity_per_sqkm:
            correction_factor:
            potential:
            min_p_max_pu:
            clip_p_max_pu:
            resource:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`snapshots_cf`, :ref:`atlite_cf`, :ref:`renewable_cf`

Inputs
------

- ``data/bundle/corine/g250_clc06_V18_5.tif``: `CORINE Land Cover (CLC) <https://land.copernicus.eu/pan-european/corine-land-cover>`_ inventory on `44 classes <https://wiki.openstreetmap.org/wiki/Corine_Land_Cover#Tagging>`_ of land use (e.g. forests, arable land, industrial, urban areas).

    .. image:: ../img/corine.png
        :scale: 33 %

- ``data/bundle/GEBCO_2014_2D.nc``: A `bathymetric <https://en.wikipedia.org/wiki/Bathymetry>`_ data set with a global terrain model for ocean and land at 15 arc-second intervals by the `General Bathymetric Chart of the Oceans (GEBCO) <https://www.gebco.net/data_and_products/gridded_bathymetry_data/>`_.

    .. image:: ../img/gebco_2019_grid_image.jpg
        :scale: 50 %

    **Source:** `GEBCO <https://www.gebco.net/data_and_products/images/gebco_2019_grid_image.jpg>`_

- ``resources/natura.tiff``: confer :ref:`natura`
- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``resources/regions_onshore.geojson``: (if not offshore wind), confer :ref:`busregions`
- ``resources/regions_offshore.geojson``: (if offshore wind), :ref:`busregions`
- ``"cutouts/" + config["renewable"][{technology}]['cutout']``: :ref:`cutout`
- ``networks/base.nc``: :ref:`base`

Outputs
-------

- ``resources/profile_{technology}.nc`` with the following structure

    ===================  ==========  =========================================================
    Field                Dimensions  Description
    ===================  ==========  =========================================================
    profile              bus, time   the per unit hourly availability factors for each node
    -------------------  ----------  ---------------------------------------------------------
    weight               bus         sum of the layout weighting for each node
    -------------------  ----------  ---------------------------------------------------------
    p_nom_max            bus         maximal installable capacity at the node (in MW)
    -------------------  ----------  ---------------------------------------------------------
    potential            y, x        layout of generator units at cutout grid cells inside the
                                     Voronoi cell (maximal installable capacity at each grid
                                     cell multiplied by capacity factor)
    -------------------  ----------  ---------------------------------------------------------
    average_distance     bus         average distance of units in the Voronoi cell to the
                                     grid node (in km)
    -------------------  ----------  ---------------------------------------------------------
    underwater_fraction  bus         fraction of the average connection distance which is
                                     under water (only for offshore)
    ===================  ==========  =========================================================

    - **profile**

    .. image:: ../img/profile_ts.png
        :scale: 33 %
        :align: center

    - **p_nom_max**

    .. image:: ../img/p_nom_max_hist.png
        :scale: 33 %
        :align: center

    - **potential**

    .. image:: ../img/potential_heatmap.png
        :scale: 33 %
        :align: center

    - **average_distance**

    .. image:: ../img/distance_hist.png
        :scale: 33 %
        :align: center

    - **underwater_fraction**

    .. image:: ../img/underwater_hist.png
        :scale: 33 %
        :align: center

Description
-----------

This script functions at two main spatial resolutions: the resolution of the
network nodes and their `Voronoi cells
<https://en.wikipedia.org/wiki/Voronoi_diagram>`_, and the resolution of the
cutout grid cells for the weather data. Typically the weather data grid is
finer than the network nodes, so we have to work out the distribution of
generators across the grid cells within each Voronoi cell. This is done by
taking account of a combination of the available land at each grid cell and the
capacity factor there.

First the script computes how much of the technology can be installed at each
cutout grid cell and each node using the `GLAES
<https://github.com/FZJ-IEK3-VSA/glaes>`_ library. This uses the CORINE land use data,
Natura2000 nature reserves and GEBCO bathymetry data.

.. image:: ../img/eligibility.png
    :scale: 50 %
    :align: center

To compute the layout of generators in each node's Voronoi cell, the
installable potential in each grid cell is multiplied with the capacity factor
at each grid cell. This is done since we assume more generators are installed
at cells with a higher capacity factor.

.. image:: ../img/offwinddc-gridcell.png
    :scale: 50 %
    :align: center

.. image:: ../img/offwindac-gridcell.png
    :scale: 50 %
    :align: center

.. image:: ../img/onwind-gridcell.png
    :scale: 50 %
    :align: center

.. image:: ../img/solar-gridcell.png
    :scale: 50 %
    :align: center

This layout is then used to compute the generation availability time series
from the weather data cutout from ``atlite``.

Two methods are available to compute the maximal installable potential for the
node (`p_nom_max`): ``simple`` and ``conservative``:

- ``simple`` adds up the installable potentials of the individual grid cells.
  If the model comes close to this limit, then the time series may slightly
  overestimate production since it is assumed the geographical distribution is
  proportional to capacity factor.

- ``conservative`` assertains the nodal limit by increasing capacities
  proportional to the layout until the limit of an individual grid cell is
  reached.

"""
import progressbar as pgb
import geopandas as gpd
import multiprocessing as mp
import xarray as xr
import pandas as pd
import numpy as np
import atlite
import logging
import rasterio as rio
import rasterio.mask
from rasterio.warp import reproject

from pypsa.geo import haversine
from shapely.geometry import LineString
from progressbar import ProgressBar
from progressbar.widgets import Percentage, SimpleProgress, Bar, Timer, ETA

from _helpers import configure_logging

logger = logging.getLogger(__name__)


def init_globals(transform_, shape_, epsg_, config_, paths_):
    """
    Import packages and defines variable for processes in multiprocess pool.
    """
    # import packages
    global rio, Resampling, np, dilation
    import rasterio as rio
    from rasterio.enums import Resampling
    import numpy as np
    from scipy.ndimage.morphology import binary_dilation as dilation
    import geopandas as gpd

    # set destination geographical data
    global dst_transform, dst_shape, dst_crs
    dst_transform = rio.Affine(*transform_)
    dst_crs = rio.crs.CRS.from_epsg(epsg_)
    dst_shape = shape_

    global config, regions
    config = config_
    paths = paths_
    regions = gpd.read_file(paths['regions']).to_crs(crs)

    # load rasters
    global gebco, clc, natura, crs, res
    natura = rasterio.open(paths['natura'])
    crs = natura.crs
    res = natura.res # reference resolution

    clc = rasterio.open(paths['corine'])
    # clc.SetProjection(gk.srs.loadSRS(3035).ExportToWkt())

    if "max_depth" in config:
        gebco = rasterio.open(paths['gebco'])
        # gebco.SetProjection(gk.srs.loadSRS(4326).ExportToWkt())





def calculate_potential(gid, save_map=None):
    """
    Calculate the potential per grid cell for one region.

    This function calculates the eligible area of the region stored in
    `path['regions']` with index `gid`. The resulting area is then projected
    onto the gridcells given in the cutout.
    """
    exclusions = []
    geom = regions.geometry.loc[[gid]]

    if config.get("natura", False):
        (natura_v,), transform = rio.mask.mask(natura, geom, crop=True, nodata=1)
        exclusions.append(natura_v)

    # nodate has value 100 (corine codes go upto 44)
    (corine_v,), corine_t = rio.mask.mask(clc, geom, crop=True, nodata=100)

    # reproject to reference (natura projection)
    kwargs = {'src_transform': corine_t, 'dst_transform': transform,
              'src_crs': crs, 'dst_crs': crs}
    corine_v, transform = reproject(corine_v, np.empty_like(natura_v), **kwargs)

    corine = config.get("corine", {})
    if isinstance(corine, list):
        corine = {'grid_codes': corine}
    if "grid_codes" in corine:
        # select codes: 1 is excluded, 0 is eligible
        ex = np.isin(corine_v, corine['grid_codes'] + [100]).astype(int)
        exclusions.append(ex)
    if corine.get("distance", 0.) > 0.:
        ex = np.isin(corine_v, corine["distance_grid_codes"] + [100]).astype(int)
        # use scipy dilation for buffer around exclusion values
        iterations = int(corine["distance"] / res[0])
        ex = dilation(ex, iterations=iterations).astype(int)
        exclusions.append(ex)


    # if "max_depth" in config:
    #     ec.excludeRasterType(gebco, (None, -config["max_depth"]))

    # # TODO compute a distance field as a raster beforehand
    # if 'max_shore_distance' in config:
    #     ec.excludeVectorType(
    #         paths["country_shapes"],
    #         buffer=config['max_shore_distance'],
    #         invert=True)
    # if 'min_shore_distance' in config:
    #     ec.excludeVectorType(
    #         paths["country_shapes"],
    #         buffer=config['min_shore_distance'])

    exclusion = (sum(exclusions) == 0).astype(float)
    return reproject(np.ones_like(exclusion), np.empty(dst_shape),
                     resampling=Resampling.bilinear,
                     src_transform=transform, dst_transform=dst_transform,
                     src_crs=crs, dst_crs=dst_crs,)[0]


if __name__ == '__main__':
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_renewable_profiles', technology='solar')
    # skip handlers due to process copying
    configure_logging(snakemake, skip_handlers=True)
    pgb.streams.wrap_stderr()
    paths = dict(snakemake.input)
    nprocesses = snakemake.config['atlite'].get('nprocesses')
    config = snakemake.config['renewable'][snakemake.wildcards.technology]
    resource = config['resource'] # pv panel config / wind turbine config
    correction_factor = config.get('correction_factor', 1.)
    capacity_per_sqkm = config['capacity_per_sqkm']
    p_nom_max_meth = config.get('potential', 'conservative')

    if correction_factor != 1.:
        logger.info(f'correction_factor is set as {correction_factor}')

    if not config.get('keep_all_available_areas', True):
        logger.warning("Argument `keep_all_available_areas` is ignored. "
                       "Continue with keeping all areas.")


    cutout = atlite.Cutout(paths['cutout'])
    minx, maxx, miny, maxy = cutout.extent
    dx = cutout.dx
    dy = cutout.dy
    # used for affine transform
    transform = [dx, 0, minx - dx / 2, 0, dy, miny - dy / 2]

    regions = gpd.read_file(paths['regions'])
    buses = pd.Index(regions.name, name='bus')

    progress = SimpleProgress(format='(%s)' %SimpleProgress.DEFAULT_FORMAT)
    widgets = [Percentage(),' ',progress,' ',Bar(),' ',Timer(),' ', ETA()]
    progressbar = ProgressBar(prefix='Compute GIS potentials: ',
                              widgets=widgets, max_value=len(regions))

    # Use the following for testing the default windows method on linux
    # mp.set_start_method('spawn')
    epsg = cutout.crs.to_epsg()
    kwargs = {'initializer': init_globals,
              'initargs': (transform, cutout.shape, epsg, config, paths),
              'maxtasksperchild': 20,
              'processes': nprocesses}
    with mp.Pool(**kwargs) as pool:
        imap = pool.imap(calculate_potential, regions.index)
        availability = list(progressbar(imap))


    coords=[('bus', buses), ('y', cutout.data.y), ('x', cutout.data.x),]
    availability = xr.DataArray(np.stack(availability), coords=coords)


    area = cutout.grid.to_crs({'proj': 'cea'}).area / 1e6
    area = xr.DataArray(area.values.reshape(cutout.shape),
                        [cutout.coords['y'], cutout.coords['x']])

    func = getattr(cutout, resource.pop('method'))
    capacity_factor = correction_factor * func(capacity_factor=True, **resource)
    capacity_potential = capacity_per_sqkm * availability.sum('bus') * area
    layout = capacity_factor * area * capacity_per_sqkm
    profile, capacities = func(matrix=availability.stack(spatial=['y','x']),
                                layout=layout, index=buses,
                                per_unit=True, return_capacity=True, **resource)

    logger.info(f"Calculating maximal capacity per bus (method '{p_nom_max_meth}')")
    if p_nom_max_meth == 'simple':
        p_nom_max = capacity_per_sqkm * availability @ area
    elif p_nom_max_meth == 'conservative':
        max_cap_factor = capacity_factor.where(availability!=0).max(['x', 'y'])
        p_nom_max = capacities / max_cap_factor
    else:
        raise AssertionError('Config key `potential` should be one of "simple" '
                        f'(default) or "conservative", not "{p_nom_max_meth}"')


    # Determine weighted average distance from substation
    layoutmatrix = (layout * availability).stack(spatial=['y','x'])
    layoutmatrix = layoutmatrix.where(capacities!=0)
    distances = haversine(regions[['x', 'y']],  cutout.grid[['x', 'y']])
    distances = layoutmatrix.copy(data=distances)
    average_distance = (layoutmatrix.weighted(distances).sum('spatial') /
                        layoutmatrix.sum('spatial'))

    ds = xr.merge([(correction_factor * profile).rename('profile'),
                    capacities.rename('weight'),
                    p_nom_max.rename('p_nom_max'),
                    capacity_potential.rename('potential'),
                    average_distance.rename('average_distance')])

    if snakemake.wildcards.technology.startswith("offwind"):
        logger.info('Calculate underwater fraction of connections.')
        offshore_shape = gpd.read_file(paths.offshore_shapes).unary_union
        underwater_fraction = []
        for i in regions.index:
            row = layoutmatrix.sel(bus=buses[i]).dropna('spatial')
            if row.data.sum() == 0:
                frac = 0
            else:
                coords = np.array([[s[1], s[0]] for s in row.spatial.data])
                centre_of_mass = coords.T @ (row.data / row.data.sum())
                line = LineString([centre_of_mass, regions.loc[i, ['x', 'y']]])
                frac = line.intersection(offshore_shape).length/line.length
            underwater_fraction.append(frac)

        ds['underwater_fraction'] = xr.DataArray(underwater_fraction, [buses])

    # select only buses with some capacity and minimal capacity factor
    ds = ds.sel(bus=((ds['profile'].mean('time') > config.get('min_p_max_pu', 0.)) &
                      (ds['p_nom_max'] > config.get('min_p_nom_max', 0.))))

    if 'clip_p_max_pu' in config:
        min_p_max_pu = config['clip_p_max_pu']
        ds['profile'] = ds['profile'].where(ds['profile'] >= min_p_max_pu, 0)

    ds.to_netcdf(snakemake.output.profile)
