..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

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

Rule ``retrieve_eurostat_data``
===============================

.. automodule:: retrieve_eurostat_data


Rule ``retrieve_jrc_idees``
===============================

.. automodule:: retrieve_jrc_idees



Rule ``retrieve_eurostat_household_data``
=========================================

.. automodule:: retrieve_eurostat_household_data


Rule ``retrieve_gas_infrastructure_data``
=========================================

.. automodule:: retrieve_gas_infrastructure_data


Rule ``retrieve_osm_data``
=========================================

.. automodule:: retrieve_osm_data

Rule ``retrieve_cutout``
============================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6382570.svg
   :target: https://doi.org/10.5281/zenodo.6382570

Cutouts are spatio-temporal subsets of the European weather data from the `ECMWF ERA5 <https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation>`__ reanalysis dataset and the `CMSAF SARAH-3 <https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>`__ solar surface radiation dataset for the year 2013, 2019 or 2023.
They have been prepared by and are for use with the `atlite <https://github.com/PyPSA/atlite>`__ tool. You can either generate them yourself using the ``build_cutouts`` rule or retrieve them directly from `zenodo <https://doi.org/10.5281/zenodo.6382570>`__ through the rule ``retrieve_cutout``.
The :ref:`tutorial` uses a smaller cutout than required for the full model (30 MB), which is also automatically downloaded.

.. note::
    To download cutouts yourself from the `ECMWF ERA5 <https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation>`__ you need to `set up the CDS API <https://cds.climate.copernicus.eu/api-how-to>`__.


**Relevant Settings**

.. code:: yaml

    tutorial:
    enable:
        build_cutout:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``cutouts/{cutout}``: weather data from either the `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`__   reanalysis weather dataset and/or `SARAH-3 <https://wui.cmsaf.eu/safira/action/viewProduktSearch>`__ satellite-based historic weather data.

.. seealso::
    For details see :mod:`build_cutout` and read the `atlite documentation <https://atlite.readthedocs.io>`__.



Rule ``retrieve_electricity_demand``
====================================

This rule downloads hourly electric load data for each country from the `OPSD platform <https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv>`__.

**Relevant Settings**

None.

**Outputs**

- ``data/electricity_demand_raw.csv``


Rule ``retrieve_cost_data``
================================

This rule downloads techno-economic assumptions from the `technology-data repository <https://github.com/pypsa/technology-data>`__.

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

Rule ``retrieve_ship_raster``
================================

This rule downloads data on global shipping traffic density from the `World Bank Data Catalogue <https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`__.

**Relevant Settings**

None.

**Outputs**

- ``data/shipdensity_global.zip``
