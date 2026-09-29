<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Data Sources

PyPSA-Eur is compiled from a variety of data sources.
The [below table](#data-inventory) provides an overview of the data sources used in PyPSA-Eur.
Different licenses apply to the data sources.

## Data Versioning {#managing_data_versions}

Many of the data sources used in PyPSA-Eur are updated regularly.
To ensure reproducibility, PyPSA-Eur uses a versioning system for data sources which
allows users to select specific versions of the data sources to use in their models.
Next to the versioning and if the license allows, most datasets are also mirrored to a
public file storage for the repository under `https://data.pypsa.org`.

!!! note
    For users, selection and control over which data sources to use is managed through the configuration file.
    See [Data Configuration](configuration.md#data_cf) for details. In most cases you just wanna stick with the latest archive
    version. Reproducibility is given even when using the `latest` tag via the
    `versions.csv`, which is version controlled.

### Understanding `versions.csv`

The file `data/versions.csv` is the central registry for all data sources and their versions.
Each row defines a specific version of a dataset with the following columns:

* `dataset`: The name of the dataset (e.g., `worldbank_urban_population`).
* `version`: The version identifier, typically following the original data source's versioning (e.g., `2025-08-14`).
* `source`: The source type - `primary` (original data source), `archive` (mirrored copy on `data.pypsa.org`), or `build` (generated from other data).
* `tags`: Space-separated tags like `latest`, `supported` or `deprecated`.
* `added`: The date when this entry was added to the registry.
* `note`: Optional notes about the dataset or version.
* `url`: The download URL for the data.

Entries to the `versions.csv` are never deleted and if a dataset was removed or is not available, the entry is marked as `deprecated`.

!!! note
    For `primary` sources, each combination of dataset and version should point to a specific version of that dataset with a unique URL.
    If the original data source does not provide versioned URLs (i.e., the URL always points to the latest data), the `version` is set to `unknown`.
    In this case, the corresponding `archive` entries do not mirror the same version but represent snapshots taken at specific points in time from that primary source.

### Updating `versions.csv`

There are two methods for updating the data versions:

1. Add a new row to `data/versions.csv`.
1. Create a new version file and reference it in the configuration under [`data→version_files`](configuration.md#data_cf).

If you are working on a fork of PyPSA-Eur, you will benefit from updating versions in a separate file as it mitigates merge conflicts when `data/versions.csv` is updated upstream.
It also clearly separates the data versions that are specific to your project from those that have been inherited from upstream.

#### Maintaining multiple version files

If you define your version updates in separate file(s), you can request them to be included at runtime by updating the `data` `version_files` key:

```yaml
data:
    version_files:
    - data/versions.csv # remove this entry if `data/my_version_overrides.csv` is a complete replacement
    - data/my_version_overrides.csv
```

New version files can be provided in CSV or YAML format.
For examples of the CSV format, see the current `data/versions.csv` file.
An equivalent YAML format example would be:

```yaml
# One list entry per equivalent CSV row
# Keys = CSV column names
- dataset: aquifer_data
    version: v1.2
    source: archive
    tags: latest supported
    added: 2025-12-02
    note: my custom note
    url: https://data.pypsa.org/workflows/eur/aquifer_data/v1.2/IHME1500_v12.zip
```

Priority is given to the definition of a dataset given later in the list.
Therefore, any datasets sharing the same `[dataset, version, source]` combination will be overwritten in the list.
Any new entries will be appended to the list.

#### Adding a new version of a dataset

=== "Updating PyPSA/PyPSA-Eur"
    If you notice that a data source has been updated and want to add the new version to PyPSA-Eur:

    1. Add a new row to `data/versions.csv` with the same `dataset` name, the new `version`, `source` set to `primary`, and the `url` pointing to the new data source.
    1. Set appropriate tags (typically `latest supported`).
    1. Update the tags of the previous version (remove `latest`, keep `supported` if still compatible).
    1. Create a pull request with your changes.
    1. Of course, any potential workflow adjustments should be considered and implemented as well.

=== "Updating your own project / soft fork"
    If you want to update a version in your own project, to one that better reflects your needs:

    1. Add a new row to your own file (e.g., `data/custom_versions.csv`) with a copy of the relevant dataset row in `data/versions.csv`.
    1. Update the `version` and `url` in this copied row, and set the `source` to `primary`.
    1. Set appropriate tags (typically `latest supported`).
    1. Add a reference to your own version file in the configuration by extending the `data→version_files` list.
    1. Of course, any potential workflow adjustments should be considered and implemented as well.

!!! note
    If the `primary` source has `version` set to `unknown` (i.e., the URL always points to the latest data) and a new version is available that has not been archived yet, please open an issue on the [PyPSA-Eur GitHub repository](https://github.com/pypsa/pypsa-eur/issues) to request an archive update.

#### Adding a new dataset

=== "Updating PyPSA/PyPSA-Eur"

    1. Add a `primary` entry to `data/versions.csv` with a new unique dataset name, version, and URL pointing to the original data source.
    1. Implement a `retrieve` rule for your dataset in `rules/retrieve.smk`.
       Take inspiration from existing rules in the file.
    1. Add the new data source to:
       * `data` section in the configuration schema `scripts/lib/validation/config/data.py`.
       * `data_inventory.csv` data inventory for PyPSA-Eur.
    1. Create a pull request with your changes.

    !!! note
        Maintainers of the repository will create the corresponding `archive` entry after reviewing your contribution.


=== "Updating your own project / soft fork"

    1. Add a `primary` entry to your own version file (e.g., `data/custom_versions.csv`) with a new unique dataset name, version, and URL pointing to the original data source.
    1. Implement a `retrieve` rule for your dataset in `rules/retrieve.smk` or your own separate retrieval rule file.
       Take inspiration from existing rules in `rules/retrieve.smk`.
    1. Add a reference to your own version file in the configuration by extending the `data→version_files` list.
    1. Add the new data source to the configuration by [extending the base schema](validation_dev.md#soft_fork_ext)

## Data inventory

{{ read_csv("doc/data_inventory.csv") }}

