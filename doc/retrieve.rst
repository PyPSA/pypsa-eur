..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _data:

###############
Retrieving Data
###############

Not all data dependencies are shipped with the git repository,
since git is not suited for handling large changing files.
Instead we provide separate data bundles which can be obtained
using the ``retrieve*`` rules.

Rule ``retrieve_databundle``
============================

.. automodule:: retrieve_databundle

Rule ``retrieve_cutout``
============================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517949.svg
   :target: https://doi.org/10.5281/zenodo.3517949

Cutouts are spatio-temporal subsets of the European weather data from the `ECMWF ERA5 <https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation>`_ reanalysis dataset and the `CMSAF SARAH-2 <https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>`_ solar surface radiation dataset for the year 2013.
They have been prepared by and are for use with the `atlite <https://github.com/PyPSA/atlite>`_ tool. You can either generate them yourself using the ``build_cutouts`` rule or retrieve them directly from `zenodo <https://doi.org/10.5281/zenodo.3517949>`__ through the rule ``retrieve_cutout``.
The :ref:`tutorial` uses a smaller cutout than required for the full model (30 MB), which is also automatically downloaded.

.. note::
    To download cutouts yourself from the `ECMWF ERA5 <https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation>`_ you need to `set up the CDS API <https://cds.climate.copernicus.eu/api-how-to>`_.


**Relevant Settings**

.. code:: yaml

    tutorial:
    enable:
        build_cutout:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``cutouts/{cutout}``: weather data from either the `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_   reanalysis weather dataset or `SARAH-2 <https://wui.cmsaf.eu/safira/action/viewProduktSearch>`_ satellite-based historic weather data.

.. seealso::
    For details see :mod:`build_cutout` and read the `atlite documentation <https://atlite.readthedocs.io>`_.


Rule ``retrieve_natura_raster``
================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4706686.svg
   :target: https://doi.org/10.5281/zenodo.4706686

This rule, as a substitute for :mod:`build_natura_raster`, downloads an already rasterized version (`natura.tiff <https://zenodo.org/record/4706686/files/natura.tiff>`_) of `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas to reduce computation times. The file is placed into the ``resources`` sub-directory.

**Relevant Settings**

.. code:: yaml

    enable:
        build_natura_raster:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``resources/natura.tiff``: Rasterized version of `Natura 2000 <https://en.wikipedia.org/wiki/Natura_2000>`_ natural protection areas to reduce computation times.

.. seealso::
    For details see :mod:`build_natura_raster`.


Rule ``retrieve_electricity_demand``
====================================

This rule downloads hourly electric load data for each country from the `OPSD platform <https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv>`_.

**Relevant Settings**

None.

**Outputs**

- ``data/load_raw.csv``


Rule ``retrieve_cost_data``
================================

This rule downloads techno-economic assumptions from the `technology-data repository <https://github.com/pypsa/technology-data>`_.

**Relevant Settings**

.. code:: yaml

    enable:
        retrieve_cost_data:

    costs:
        year:
        version:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`costs_cf`

**Outputs**

- ``resources/costs.csv``

Rule ``retrieve_irena``
================================

.. automodule:: retrieve_irena

Rule ``retrieve_ship_raster``
================================

This rule downloads data on global shipping traffic density from the `World Bank Data Catalogue <https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`_.

**Relevant Settings**

None.

**Outputs**

- ``data/shipdensity_global.zip``


Rule ``retrieve_sector_databundle``
====================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5546516.svg
   :target: https://doi.org/10.5281/zenodo.5546516

In addition to the databundle required for electricity-only studies,
another databundle is required for modelling sector-coupled systems.
The size of this data bundle is around 640 MB.
