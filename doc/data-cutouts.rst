.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _cutouts:

############
Weather Data
############

Cutouts are spatio-temporal subsets of the European weather data from the `ECMWF ERA5 <https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation>`__ reanalysis dataset and the `CMSAF SARAH-3 <https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>`__ solar surface radiation dataset.
They have been prepared by and are for use with the `atlite <https://github.com/PyPSA/atlite>`__ tool.
The :ref:`tutorial` uses a smaller cutout than required for the full model (30 MB), which is also automatically downloaded.

There are two ways to obtain cutouts:

- **Retrieve pre-built cutouts** from the archive using the ``retrieve_cutout`` rule (default). This downloads ready-made cutouts from ``data.pypsa.org``.
- **Build cutouts from scratch** using the ``build_cutout`` rule. This requires access to the `CDS API <https://cds.climate.copernicus.eu/api-how-to>`__ to download ERA5 data directly.

.. seealso::
    For building your own cutouts, see :mod:`build_cutout` and the `atlite documentation <https://atlite.readthedocs.io>`__.

The following pre-built cutouts are available for download under
``https://data.pypsa.org/workflows/cutout/<version>/<cutout>.nc``
(click on a cutout name below to download directly).

.. _available-cutouts-start:

**v1.0**

.. list-table::
   :header-rows: 1
   :widths: 70 30

   * - Cutout
     - Size
   * - `be-03-2013-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/be-03-2013-era5.nc>`__
     - 11.1 MB
   * - `dach-03-2013-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/dach-03-2013-sarah3-era5.nc>`__
     - 35.1 MB
   * - `europe-1995-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-1995-sarah3-era5.nc>`__
     - 6.2 GB
   * - `europe-1996-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-1996-sarah3-era5.nc>`__
     - 6.5 GB
   * - `europe-2008-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2008-sarah3-era5.nc>`__
     - 6.2 GB
   * - `europe-2009-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2009-sarah3-era5.nc>`__
     - 6.2 GB
   * - `europe-2010-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2010-sarah3-era5.nc>`__
     - 6.1 GB
   * - `europe-2012-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2012-sarah3-era5.nc>`__
     - 6.6 GB
   * - `europe-2013-03-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2013-03-sarah3-era5.nc>`__
     - 140.5 MB
   * - `europe-2013-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2013-sarah3-era5.nc>`__
     - 6.1 GB
   * - `europe-2019-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2019-sarah3-era5.nc>`__
     - 6.1 GB
   * - `europe-2020-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2020-sarah3-era5.nc>`__
     - 6.6 GB
   * - `europe-2021-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2021-sarah3-era5.nc>`__
     - 6.1 GB
   * - `europe-2023-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2023-sarah3-era5.nc>`__
     - 6.1 GB
   * - `europe-2024-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2024-sarah3-era5.nc>`__
     - 6.7 GB
   * - `europe-2025-sarah3-era5.nc <https://data.pypsa.org/workflows/cutout/v1.0/europe-2025-sarah3-era5.nc>`__
     - 6.7 GB
.. _available-cutouts-end:

**Relevant Settings**

.. code:: yaml

    atlite:
      default_cutout:
      cutouts:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`atlite_cf` and :ref:`data_cf`.

**Outputs**

- ``cutouts/{cutout}``: weather data from either the `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`__   reanalysis weather dataset and/or `SARAH-3 <https://wui.cmsaf.eu/safira/action/viewProduktSearch>`__ satellite-based historic weather data.


