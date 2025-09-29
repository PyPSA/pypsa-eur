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

######################
Managing Data Versions
######################

Many of the data sources used in PyPSA-Eur are updated regularly.
To ensure reproducibility, PyPSA-Eur uses a versioning system for data sources which allows users to select specific versions of the data sources to use in their models.

.. note::
    For users, selection and control over which is managed through the configuration file.
    See :ref:`_data_cf` for details.

##########################################
Creating a new version of the data sources
##########################################

To create a new version of the data sources, you can use the helper script in ``scripts/create_zenodo_deposition_cli.py``.
Here are the steps that this script helps you to navigate:

1. Locate the data for the new version and place it under ``data/<dataset_name>/archive/<version>/``.
   E.g. for creating a new version ``2029-01-01`` of the ``worldbank_population`` dataset, place the data into a folder named ``data/worldbank_population/archive/2029-01-01/``.
   We follow the versioning names of the original dataset, so make sure to use the same version name as the original dataset.
2. If you want to use the script, run it now. It will guide you through the process outlined below.
3. Create a new Zenodo deposition for the new version of the data. 
   You can do this by visiting the Zenodo website and creating a new deposition.
   * All relevant metadata, such as the title, description, and keywords based on the previous version of the dataset and make changes as needed.
   * Make sure to add the previous version Zenodo deposit as a related identifier with the relation "is new version of", "doi" as the identifier scheme, and the DOI of the previous version as the identifier.
   * Upload the data files from the ``data/<dataset_name>/archive/<version>/`` folder to the Zenodo deposition.
4. Once the deposition is complete, publish it on Zenodo.
5. Update ``data/versions.csv`` to include the new version of the dataset.
   * Create a new row in the CSV file based on the previous version, updating the version number and Zenodo URL
   * Make sure to tag this new version with the tags ``['latest', 'supported']``.
   * Remove the ``latest`` tag from the previous version.
   * If the previous version is no longer supported or outdated, remove the ``supported`` tag and add the ``deprecated`` tag.
6. Commit the changes to the repository and create a pull request.

########################
Adding a new data source
########################

The process of adding a new data source is similar to creating a new version of an existing data source, with some additional steps.
It is also possible to use the helper script in ``scripts/create_zenodo_deposition_cli.py`` to guide you through the process.

1. Create a new folder for the new data source in the ``data/`` directory, e.g. ``data/my_new_data_source/``.
2. Place the data files for the new data source in the ``data/<dataset_name>/archive/<version>/`` folder.
   We follow the versioning names of the original dataset, so make sure to use the same version name as the original dataset.
3. If you want to use the script, run it now. It will guide you through the process outlined below.
4. Create a new Zenodo deposition for the new data source.
   You can do this by visiting the Zenodo website and creating a new deposition.
   * Add all relevant metadata, such as the title, authors, description, keywords.
   * When adding the license, make sure that the license of the dataset is compatible with redistribution (i.e. uploading to Zenodo). Most of our data is originally CC-BY-4.0 licensed. If you have doubts about the license, reach out to the maintainers.
   * Make sure to set the version name to be the same as the version name used for the data files, e.g. ``2029-01-01``.
   * Upload the data files from the ``data/<dataset_name>/archive/<version>/`` folder to the Zenodo deposition. 
5. Once the deposition is complete, publish it on Zenodo.
6. Update ``data/versions.csv`` to include the new data source.
   * Create a new row in the CSV file with the following columns:
     * ``dataset``: The name of the dataset as used in the folder name, e.g. ``my_new_data_source``.
     * ``source``: The source of the dataset. For Zenodo uploads the source is by definition ``archive``.
     * ``version``: The version name of the dataset as used in the folder name, e.g. ``2029-01-01``.
     * ``tags``: A list of tags for the dataset. Make sure to include ``latest`` and ``supported`` tags.
     * ``url``: The link to the Zenodo deposition of the dataset, e.g. ``https://zenodo.org/record/<zenodo_id>``. Check whether the respective ``retrieve_<dataset_name>`` rule in ``rules/retrieve.smk`` requires a direct download link or the link to the Zenodo record.
     * ``note``: An optional note about the dataset.
7. Implement a ``retrieve`` rule for your dataset in ``rules/retrieve.smk``.
   This rule should download the data from the Zenodo deposition and place it in the ``data/<dataset_name>/archive/<version>/`` folder.
   Take inspiration from existing rules in the file, e.g. the ``rule retrieve_worldbank_urban_population``.
8. Add an additional rule for the ``primary`` source of the data, i.e. the original source of the data.
   You may be able to use the same rule for ``archive`` and ``primary`` sources, but sometimes dedicated rules are needed.
   Create an entry in ``data/versions.csv`` for the ``primary`` source as well, with the URL pointing to the original source of the data.
   Again, take inspiration from existing rules in the file, e.g. the ``rule retrieve_worldbank_urban_population_primary``.
   This rule will also help us in the future to update to new versions of the data set.
9. Add the new data source to the
    * ``data`` section in the configuration file ``config/config.default.yaml``
    * ``doc/configtables/data.csv`` for the documentation
    * ``data_sources.rst`` data inventory for PyPSA-Eur

==============
Data inventory
==============

.. csv-table::
   :header-rows: 1
   :class: longtable
   :widths: auto
   :file: data_inventory.csv