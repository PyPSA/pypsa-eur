# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from scripts.lib.validation.config._base import ConfigModel
from scripts.lib.validation.config.adjustments import AdjustmentsConfig
from scripts.lib.validation.config.atlite import AtliteConfig
from scripts.lib.validation.config.biomass import BiomassConfig
from scripts.lib.validation.config.clustering import ClusteringConfig
from scripts.lib.validation.config.co2_budget import Co2BudgetConfig
from scripts.lib.validation.config.conventional import ConventionalConfig
from scripts.lib.validation.config.costs import CostsConfig
from scripts.lib.validation.config.countries import CountriesConfig
from scripts.lib.validation.config.data import DataConfig
from scripts.lib.validation.config.electricity import ElectricityConfig
from scripts.lib.validation.config.enable import EnableConfig
from scripts.lib.validation.config.energy import EnergyConfig
from scripts.lib.validation.config.existing_capacities import ExistingCapacitiesConfig
from scripts.lib.validation.config.foresight import ForesightConfig
from scripts.lib.validation.config.industry import IndustryConfig
from scripts.lib.validation.config.lines import LinesConfig
from scripts.lib.validation.config.links import LinksConfig
from scripts.lib.validation.config.load import LoadConfig
from scripts.lib.validation.config.overpass_api import OverpassApiConfig
from scripts.lib.validation.config.pypsa_eur import PypsaEurConfig
from scripts.lib.validation.config.renewable import RenewableConfig
from scripts.lib.validation.config.run import RunConfig
from scripts.lib.validation.config.scenario import ScenarioConfig
from scripts.lib.validation.config.sector import SectorConfig
from scripts.lib.validation.config.snapshots import SnapshotsConfig
from scripts.lib.validation.config.solar_thermal import SolarThermalConfig
from scripts.lib.validation.config.solving import SolvingConfig
from scripts.lib.validation.config.transformers import TransformersConfig
from scripts.lib.validation.config.transmission_projects import (
    TransmissionProjectsConfig,
)


class LoggingConfig(ConfigModel):
    """Configuration for top level `logging` settings."""

    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = Field(
        "INFO",
        description="Restrict console outputs to all infos, warning or errors only",
    )
    format: str = Field(
        "%(levelname)s:%(name)s:%(message)s",
        description="Custom format for log messages. See `LogRecord <https://docs.python.org/3/library/logging.html#logging.LogRecord>`_ attributes.",
    )


class RemoteConfig(ConfigModel):
    """Configuration for top level `remote` settings."""

    ssh: str = Field(
        "",
        description="Optionally specify the SSH of a remote cluster to be synchronized.",
    )
    path: str = Field(
        "",
        description="Optionally specify the file path within the remote cluster to be synchronized.",
    )


class ConfigSchema(BaseModel):
    """
    Combined configuration schema for PyPSA-EUR.
    """

    # TODO Change to extra='forbid' once schema covers all config options
    # For soft-forks it is recommended to either extend the schema for full config
    # coverage or allow extra fields with extra='allow'
    model_config = ConfigDict(extra="allow", title="PyPSA-Eur Configuration")

    _name: str = "default"
    """internal attribute to track the config filename following the application of config updates"""

    # Top-level fields (from TopLevelConfig)
    version: str = Field(
        "v2025.07.0", description="Version of PyPSA-Eur. Descriptive only."
    )
    tutorial: bool = Field(
        False,
        description="Switch to retrieve the tutorial data set instead of the full data set.",
    )
    logging: LoggingConfig = Field(
        default_factory=LoggingConfig,
        description="Logging configuration for the workflow",
    )
    remote: RemoteConfig = Field(
        default_factory=RemoteConfig,
        description="Configuration for remote workflow execution",
    )

    run: RunConfig = Field(
        default_factory=RunConfig,
        description="Run configuration for PyPSA-EUR workflow execution.",
    )
    foresight: ForesightConfig = Field(
        default_factory=ForesightConfig,
        description="Foresight mode for the optimization. See Foresight Options for detailed explanations.",
    )
    scenario: ScenarioConfig = Field(
        default_factory=ScenarioConfig,
        description="Scenario configuration defining wildcards for the workflow.",
    )
    countries: CountriesConfig = Field(
        default_factory=CountriesConfig,
        description="European countries defined by their Two-letter country codes (ISO 3166-1) which should be included in the energy system model.",
    )
    snapshots: SnapshotsConfig = Field(
        default_factory=SnapshotsConfig,
        description="Configuration for the time period snapshots.",
    )
    enable: EnableConfig = Field(
        default_factory=EnableConfig,
        description="Flags to enable/disable workflow features.",
    )
    co2_budget: Co2BudgetConfig | None = Field(
        default_factory=Co2BudgetConfig,
        description="CO2 budget as fraction of 1990 emissions per planning horizon year.",
    )
    electricity: ElectricityConfig = Field(
        default_factory=ElectricityConfig,
        description="Electricity sector configuration.",
    )
    atlite: AtliteConfig = Field(
        default_factory=AtliteConfig,
        description="Atlite cutout configuration for weather data.",
    )
    renewable: RenewableConfig = Field(
        default_factory=RenewableConfig,
        description="Renewable energy technologies configuration.",
    )
    conventional: ConventionalConfig = Field(
        default_factory=ConventionalConfig,
        description="Conventional power plants configuration.",
    )
    lines: LinesConfig = Field(
        default_factory=LinesConfig,
        description="Transmission lines configuration.",
    )
    links: LinksConfig = Field(
        default_factory=LinksConfig,
        description="HVDC links configuration.",
    )
    transmission_projects: TransmissionProjectsConfig = Field(
        default_factory=TransmissionProjectsConfig,
        description="Transmission projects configuration.",
    )
    transformers: TransformersConfig = Field(
        default_factory=TransformersConfig,
        description="Transformers configuration.",
    )
    load: LoadConfig = Field(
        default_factory=LoadConfig,
        description="Electrical load configuration.",
    )
    pypsa_eur: PypsaEurConfig = Field(
        default_factory=PypsaEurConfig,
        description="PyPSA-Eur component filtering configuration.",
    )
    energy: EnergyConfig = Field(
        default_factory=EnergyConfig,
        description="Energy totals configuration.",
    )
    biomass: BiomassConfig = Field(
        default_factory=BiomassConfig,
        description="Biomass configuration.",
    )
    solar_thermal: SolarThermalConfig = Field(
        default_factory=SolarThermalConfig,
        description="Solar thermal configuration.",
    )
    existing_capacities: ExistingCapacitiesConfig = Field(
        default_factory=ExistingCapacitiesConfig,
        description="Existing capacities grouping configuration.",
    )
    sector: SectorConfig = Field(
        default_factory=SectorConfig,
        description="Sector coupling configuration.",
    )
    industry: IndustryConfig = Field(
        default_factory=IndustryConfig,
        description="Industry sector configuration.",
    )
    costs: CostsConfig = Field(
        default_factory=CostsConfig,
        description="Cost assumptions configuration.",
    )
    clustering: ClusteringConfig = Field(
        default_factory=ClusteringConfig,
        description="Network clustering configuration.",
    )
    adjustments: AdjustmentsConfig = Field(
        default_factory=AdjustmentsConfig,
        description="Network adjustments configuration.",
    )
    solving: SolvingConfig = Field(
        default_factory=SolvingConfig,
        description="Solver and optimization configuration.",
    )
    data: DataConfig = Field(
        default_factory=DataConfig,
        description="Data source configuration.",
    )
    overpass_api: OverpassApiConfig = Field(
        default_factory=OverpassApiConfig,
        description="Overpass API configuration for OSM data retrieval.",
    )

    @model_validator(mode="before")
    @classmethod
    def check_no_secrets_section(cls, data):
        """Prevent secrets from being stored in config."""
        if isinstance(data, dict) and "secrets" in data:
            raise ValueError(
                "The 'secrets:' section is no longer supported in config to avoid "
                "leaking credentials. Use environment variables instead (e.g., "
                "CORINE_API_TOKEN). You can set these in a .env file in the project root."
            )
        return data
