# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Data source configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#data
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


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
    bidding_zones_electricitymaps: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Electricitymaps bidding zones data source configuration.",
    )
    bidding_zones_entsoepy: _DataSourceConfig = Field(
        default_factory=_DataSourceConfig,
        description="Entsoepy bidding zones data source configuration.",
    )
