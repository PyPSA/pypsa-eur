.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _data:

###############
Retrieving Data
###############

Not all data dependencies are shipped with the git repository, since git is not suited for handling large changing files.
Instead we use separate steps in the workflow (``rules`` executed by ``snakemake``) to download external data using the ``retrieve_<dataset>`` rules.

Data is generally retrieved in a version-controlled manner, enabling control over input data versions, reproducibility and consistency of modelling runs.
The rules download data into subfolders in the `data/` directory, following the structure 
``data/{dataset}/{source}/{version}``, e.g. ``data/jrc_idees/primary/March-2025-V1/``.
Which specific data version is retrieve can be controlled in the `data configuration <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#data>`__ .

Below some specific ``retrieve_<dataset>`` rules are documented.
For more information on the datasets retrieved, see the `data sources <https://pypsa-eur.readthedocs.io/en/latest/data_sources.html>`__ and *Data inventory* section there in the documentation.

Rule ``retrieve_bidding_zones``
=========================================

.. automodule:: retrieve_bidding_zones

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
    :ref:`tutorial_cf` and :ref:`enable_cf`.

**Outputs**

- ``cutouts/{cutout}``: weather data from either the `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`__   reanalysis weather dataset and/or `SARAH-3 <https://wui.cmsaf.eu/safira/action/viewProduktSearch>`__ satellite-based historic weather data.

.. seealso::
    For details see :mod:`build_cutout` and read the `atlite documentation <https://atlite.readthedocs.io>`__.


Rule ``retrieve_cost_data``
================================

This rule downloads techno-economic assumptions from the `technology-data repository <https://github.com/pypsa/technology-data>`__.

**Relevant Settings**

.. code:: yaml

    costs:
        year:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`costs_cf`

**Outputs**

- ``data/costs/primary/{version}/costs_{year}.csv``