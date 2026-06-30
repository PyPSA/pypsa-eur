# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Data source configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#data
"""

from datetime import date
from pathlib import Path
from typing import Literal

import pandas as pd
from natsort import natsort_keygen
from pandera.pandas import Check, Column, DataFrameSchema
from pydantic import BaseModel, Field, FilePath, field_validator

from scripts.lib.validation.config._base import ConfigModel

VALID_SOURCES = ["primary", "archive", "build"]  # Order defines sort priority

VALID_TAGS = {
    "latest",
    "supported",
    "not-supported",
    "deprecated",
    "might-work",
    "not-tested",
    "broken-link",
}

not_empty = [Check.str_length(min_value=1), Check.str_matches(r"\S")]
valid_tags = Check(lambda s: all(t in VALID_TAGS for t in s.split()), element_wise=True)
url_safe = Check.str_matches(
    r"^[a-z0-9_\-\.]+$",
    error="Version must be URL-safe (only alphanumeric, hyphen, underscore, dot)",
)


def sort_versions(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by dataset (asc), version (desc, natural), source (as predefined)."""
    natsort_key = natsort_keygen(key=str.casefold)
    df = df.copy()
    df["source"] = pd.Categorical(df["source"], categories=VALID_SOURCES, ordered=True)
    df = df.sort_values(
        by=["dataset", "version", "source"],
        key=lambda c: c.map(natsort_key) if c.name != "source" else c,
        ascending=[True, False, True],
    ).reset_index(drop=True)
    df["source"] = df["source"].astype(str)
    df = df[VersionsSchema.columns.keys()]  # Ensure column order
    return df


is_sorted = Check(lambda df: df.equals(sort_versions(df)), error="Data must be sorted")
archive_has_url = Check(
    lambda df: df.loc[df["source"] == "archive", "url"].str.len().gt(0).all(),
    error="Archive entries must have a URL",
)
one_latest_per_dataset_source = Check(
    lambda df: (
        df[df["tags"].str.contains("latest")]
        .groupby(["dataset", "source"])
        .size()
        .eq(1)
        .all()
    ),
    error="Exactly one 'latest' tag required per dataset/source combination",
)
latest_same_version_across_sources = Check(
    lambda df: (
        df[(df["tags"].str.contains("latest")) & (df["version"] != "unknown")]
        .groupby("dataset")["version"]
        .nunique()
        .le(1)
        .all()
    ),
    error="All 'latest' entries for a dataset must have the same version across sources (excluding 'unknown')",
)
VersionsSchema = DataFrameSchema(
    {
        "dataset": Column(str, not_empty, nullable=False),
        "version": Column(str, not_empty + [url_safe], nullable=False),
        "source": Column(str, Check.isin(VALID_SOURCES)),
        "tags": Column(str, valid_tags, nullable=False),
        "added": Column(
            str,
            Check.str_matches(r"^\d{4}-\d{2}-\d{2}$"),
            nullable=False,
            coerce=True,
            default=date.today().isoformat(),
        ),
        "note": Column(str, nullable=True),
        "url": Column(
            str,
            Check.str_matches(
                r'^(https?://[^:\s<>"|*]*)?$',
                error='URL must start with http(s):// and not contain colons (use %3A), spaces, or Windows-invalid characters (<>"|*).',
            ),
            nullable=True,
        ),
    },
    checks=[
        is_sorted,
        archive_has_url,
        one_latest_per_dataset_source,
        latest_same_version_across_sources,
    ],
    coerce=True,
    strict=True,
    ordered=True,
    unique=["dataset", "version", "source"],
)


class _DataSourceConfig(ConfigModel):
    """Configuration for a single data source."""

    source: Literal["archive", "primary", "build"] = Field(
        "archive",
        description="Source of the data. 'archive' retrieves pre-built data, 'primary' retrieves from primary source.",
    )
    version: str = Field(
        "latest",
        description="Version of the data to use. Uses the specific 'version' for the selected 'source' or the dataset tagged 'latest' for this source.",
    )


class DataConfig(BaseModel):
    """Configuration for `data` settings."""

    @field_validator("version_files")
    @classmethod
    def check_version_files_are_correct_suffix(
        cls, v: list[FilePath]
    ) -> list[FilePath]:
        for path in v:
            if path.suffix.lower() not in {".csv", ".yaml", ".yml"}:
                raise ValueError(f"Version file '{path}' must be a CSV or YAML file.")
        return v

    version_files: list[FilePath] = Field(
        default_factory=lambda: [Path("data/versions.csv")],
        description="""
        List of paths to version files.
        If multiple paths are provided, they will be merged with priority given to the later paths in the list.
        This allows for overriding default versions with custom versions, without modifying the default version file.

        All filepaths must be relative to the project root or as an absolute path and point to one of:
        - a CSV file.
        - a YAML file with a list of version entries, which will be converted internally to CSV format.
          This allows for easier editing of versions with comments and multi-line notes.
          This is equivalent to calling `df.to_dict(orient='records')` on a dataframe with the same structure as the CSV version files, and then dumping to YAML.

        The CSV columns / YAML fields must be:
        - `dataset`: The name of the dataset (e.g., "hotmaps_industrial_sites", "enspreso_biomass", etc.).
        - `version`: The version of the dataset to use (e.g., "latest", "v1.0", etc.).
        - `source`: The source of the dataset (e.g., "archive", "primary", "build").
          "archive" retrieves data from an organisation data bucket (e.g. `data.pypsa.org`).
          "primary" retrieves data from the primary source (e.g. Eurostat, OSM, etc.).
          "build" retrieves data from the result of a build script within the PyPSA-Eur workflow itself, rather than a remote source.
        - `tags`: space separated tags for the dataset, which will become boolean columns when loaded into the workflow.
          "supported" must be a tag on a version if it is supported by the workflow.
          "latest" must be a tag on the version of each `source` that should be used when `version` is set to "latest".
        - `added`: The date when the dataset version was added (e.g., "2024-01-01").
        - `note`: [Optional] notes about the dataset version.
        - `url`: URL to the dataset version. Optional if data `source` is "build", otherwise required.
        """,
    )
    hotmaps_industrial_sites: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Hotmaps industrial sites data source configuration.",
    )
    enspreso_biomass: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Enspreso biomass data source configuration.",
    )
    osm: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="OSM data source configuration.",
    )
    worldbank_urban_population: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="World Bank urban population data source configuration.",
    )
    worldbank_commodity_prices: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="World Bank commodity prices data source configuration.",
    )
    gem_europe_gas_tracker: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="GEM Europe Gas Tracker data source configuration.",
    )
    gem_gcct: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="GEM Global Cement and Concrete Tracker data source configuration.",
    )
    instrat_co2_prices: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="Instrat CO2 prices data source configuration.",
    )
    co2stop: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="CO2Stop data source configuration.",
    )
    nitrogen_statistics: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Nitrogen statistics data source configuration.",
    )
    eu_nuts2013: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="EU NUTS 2013 data source configuration.",
    )
    eu_nuts2021: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="EU NUTS 2021 data source configuration.",
    )
    eurostat_balances: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Eurostat balances data source configuration.",
    )
    eurostat_household_balances: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Eurostat household balances data source configuration.",
    )
    wdpa: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="WDPA data source configuration.",
    )
    wdpa_marine: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="WDPA Marine data source configuration.",
    )
    luisa_land_cover: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="LUISA land cover data source configuration.",
    )
    jrc_idees: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="JRC IDEES data source configuration.",
    )
    scigrid_gas: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="SciGRID Gas data source configuration.",
    )
    seawater_temperature: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Seawater temperature data source configuration.",
    )
    swiss_energy_balances: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="Swiss energy balances data source configuration.",
    )
    synthetic_electricity_demand: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="Synthetic electricity demand data source configuration.",
    )
    opsd_electricity_demand: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="archive"),
        description="OPSD electricity demand data source configuration.",
    )
    entsoe_electricity_demand: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="archive"),
        description="ENTSO-E electricity demand data source configuration.",
    )
    neso_electricity_demand: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="archive"),
        description="NESO electricity demand data source configuration.",
    )
    copernicus_land_cover: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="Copernicus land cover data source configuration.",
    )
    ship_raster: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Ship raster data source configuration.",
    )
    eez: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="EEZ data source configuration.",
    )
    nuts3_population: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="NUTS3 population data source configuration.",
    )
    gdp_per_capita: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="GDP per capita data source configuration.",
    )
    population_count: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Population count data source configuration.",
    )
    ghg_emissions: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="GHG emissions data source configuration.",
    )
    gebco: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="GEBCO data source configuration.",
    )
    attributed_ports: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Attributed ports data source configuration.",
    )
    corine: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="CORINE data source configuration.",
    )
    emobility: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="E-mobility data source configuration.",
    )
    h2_salt_caverns: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="H2 salt caverns data source configuration.",
    )
    lau_regions: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="LAU regions data source configuration.",
    )
    aquifer_data: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Aquifer data source configuration.",
    )
    osm_boundaries: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="OSM boundaries data source configuration.",
    )
    gem_gspt: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="GEM GSPT data source configuration.",
    )
    tyndp: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="TYNDP data source configuration.",
    )
    powerplants: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Powerplants data source configuration.",
    )
    costs: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Costs data source configuration.",
    )
    country_runoff: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Country runoff data source configuration.",
    )
    country_hdd: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Country HDD data source configuration.",
    )
    natura: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Natura data source configuration.",
    )
    bfs_road_vehicle_stock: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="BFS road vehicle stock data source configuration.",
    )
    bfs_gdp_and_population: _DataSourceConfig = Field(
        default_factory=lambda: _DataSourceConfig(source="primary"),
        description="BFS GDP and population data source configuration.",
    )
    mobility_profiles: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Mobility profiles data source configuration.",
    )
    cutout: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Cutout data source configuration.",
    )
    dh_areas: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="District heating areas data source configuration.",
    )
    geothermal_heat_utilisation_potentials: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Geothermal heat utilisation potentials data source configuration.",
    )
    jrc_ardeco: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="JRC ARDECO data source configuration.",
    )
    jrc_energy_atlas: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="JRC Energy Atlas data source configuration.",
    )
    desnz_electricity_consumption: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="DESNZ (UK Department for Energy Security and Net Zero) electricity consumption data source configuration.",
    )
    ons_lad: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="ONS (Office for National Statistics) Local Authority District data source configuration.",
    )
    bidding_zones_electricitymaps: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Electricitymaps bidding zones data source configuration.",
    )
    bidding_zones_entsoepy: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Entsoepy bidding zones data source configuration.",
    )
