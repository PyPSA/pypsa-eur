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

See :ref:`cutouts`.


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