.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

#############
Data Sources
#############

PyPSA-Eur is compiled from a variety of data sources. The following table provides an
overview of the data sources used in PyPSA-Eur. Different licenses apply to the
data sources.

.. toctree::
   :maxdepth: 1

   ../data-base-network
   ../data-repos

.. _managing_data_versions:

#################
Data Versioning
#################

Many of the data sources used in PyPSA-Eur are updated regularly.
To ensure reproducibility, PyPSA-Eur uses a versioning system for data sources which
allows users to select specific versions of the data sources to use in their models.
Next to the versioning and if the license allows, most datasets are also mirrored to a
public file storage for the repository under ``https://data.pypsa.org``.

.. note::
    For users, selection and control over which data sources to use is managed through the configuration file.
    See :ref:`data_cf` for details. In most cases you just wanna stick with the latest archive
    version. Reproducibility is given even when using the ``latest`` tag via the
    ``versions.csv``, which is version controlled.

*****************************
Understanding ``versions.csv``
*****************************

The file ``data/versions.csv`` is the central registry for all data sources and their versions.
Each row defines a specific version of a dataset with the following columns:

* ``dataset``: The name of the dataset (e.g., ``worldbank_urban_population``).
* ``version``: The version identifier, typically following the original data source's versioning (e.g., ``2025-08-14``).
* ``source``: The source type - ``primary`` (original data source), ``archive`` (mirrored copy on ``data.pypsa.org``), or ``build`` (generated from other data).
* ``tags``: Space-separated tags like ``latest``, ``supported`` or ``deprecated``.
* ``added``: The date when this entry was added to the registry.
* ``note``: Optional notes about the dataset or version.
* ``url``: The download URL for the data.

Entries to the ``versions.csv`` are never deleted and if a dataset was removed or is not available, the entry is marked as ``deprecated``.

.. note::
    For ``primary`` sources, each combination of dataset and version should point to a specific version of that dataset with a unique URL.
    If the original data source does not provide versioned URLs (i.e., the URL always points to the latest data), the ``version`` is set to ``unknown``.
    In this case, the corresponding ``archive`` entries do not mirror the same version but represent snapshots taken at specific points in time from that primary source.

*******************************
Adding a new version of a dataset
*******************************

If you notice that a data source has been updated and want to add the new version to PyPSA-Eur:

1. Add a new row to ``data/versions.csv`` with the same ``dataset`` name, the new ``version``, ``source`` set to ``primary``, and the ``url`` pointing to the original data source.
2. Set appropriate tags (typically ``latest supported``).
3. Update the tags of the previous version (remove ``latest``, keep ``supported`` if still compatible).
4. Create a pull request with your changes.
5. Of course, any potential workflow adjustments should be considered and implemented as well.

.. note::
    If the ``primary`` source has ``version`` set to ``unknown`` (i.e., the URL always points to the latest data) and a new version is available that has not been archived yet, please open an issue on the `PyPSA-Eur GitHub repository <https://github.com/pypsa/pypsa-eur/issues>`_ to request an archive update.

*********************
Adding a new dataset
*********************

To add a completely new data source to PyPSA-Eur:

1. Add a ``primary`` entry to ``data/versions.csv`` with a new unique dataset name, version, and URL pointing to the original data source.
2. Implement a ``retrieve`` rule for your dataset in ``rules/retrieve.smk``.
   Take inspiration from existing rules in the file.
3. Add the new data source to:

   * ``data`` section in the pydantic schema ``scripts/lib/validation/config/data.py``
   * ``data_inventory.csv`` data inventory for PyPSA-Eur

4. Create a pull request with your changes.

.. note::
    Maintainers of the repository will create the corresponding ``archive`` entry after reviewing your contribution.

##############
Data inventory
##############

.. csv-table::
   :header-rows: 1
   :class: longtable
   :widths: auto
   :file: data_inventory.csv
