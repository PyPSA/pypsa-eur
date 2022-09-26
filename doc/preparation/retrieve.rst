..
  SPDX-FileCopyrightText: 2019-2022 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _data:

Rules ``retrieve*``
=============================

Not all data dependencies are shipped with the git repository,
since git is not suited for handling large changing files.
Instead we provide separate data bundles which can be obtained
using the ``retrieve*`` rules.

Rule ``retrieve_databundle``
----------------------------

.. automodule:: retrieve_databundle

Rule ``retrieve_cutout``
------------------------

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
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``cutouts/{cutout}``: weather data from either the `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_   reanalysis weather dataset or `SARAH-2 <https://wui.cmsaf.eu/safira/action/viewProduktSearch>`_ satellite-based historic weather data.

.. seealso::
    For details see :mod:`build_cutout` and read the `atlite documentation <https://atlite.readthedocs.io>`_.
