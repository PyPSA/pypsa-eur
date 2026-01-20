# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Config validation for PyPSA-EUR.

The schema is exported to both `config/config.default.yaml` and `config/schema.json`.
The json schema is also contributed to the schemastore.org and matches
`**/pypsa-eur*/config/*.yaml` to get IDE support without additional configuration.
"""

import re
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator
from ruamel.yaml import YAML

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
    co2_budget: Co2BudgetConfig = Field(
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


def validate_config(config: dict) -> ConfigSchema:
    """Validate config dict against schema."""
    return ConfigSchema(**config)


def generate_config_defaults(path: str = "config/config.default.yaml") -> dict:
    """Generate config defaults YAML file and return the defaults dict."""
    from ruamel.yaml.comments import CommentedMap

    def convert_to_field_name(key: str) -> str:
        """Convert dash-case to snake_case for field lookup."""
        return key.replace("-", "_")

    # by_alias is needed to export dash-case instead of snake_case (which are some set aliases)
    # the goal should be to use snake_case consistently
    defaults = ConfigSchema().model_dump(by_alias=True)

    # Create YAML instance with custom settings
    yaml_writer = YAML()
    yaml_writer.version = (1, 1)  # Make sure to quote boolean-looking strings
    yaml_writer.default_flow_style = False
    yaml_writer.width = 4096  # Avoid line wrapping
    yaml_writer.indent(mapping=2, sequence=2, offset=0)

    # Custom string representer for controlling quote style
    def str_representer(dumper, data):
        """Use block style for multiline, quotes for special chars, plain otherwise."""
        TAG = "tag:yaml.org,2002:str"
        if "\n" in data:
            return dumper.represent_scalar(TAG, data, style="|")
        if data == "" or any(c in data for c in ":{}[]&*#?|-<>=!%@"):
            return dumper.represent_scalar(TAG, data, style='"')
        return dumper.represent_scalar(TAG, data, style="")

    yaml_writer.representer.add_representer(str, str_representer)

    # Create a CommentedMap to add comments
    data = CommentedMap()

    # Add yaml-language-server comment at the very top (before first key)
    data.yaml_set_start_comment("yaml-language-server: $schema=./schema.json")

    for key, value in defaults.items():
        data[key] = value

        field_name = convert_to_field_name(key)
        docs_url = f"https://pypsa-eur.readthedocs.io/en/latest/configuration.html#{field_name}"
        data.yaml_set_comment_before_after_key(key, before=f"\ndocs in {docs_url}")

    # Write to file
    with open(path, "w") as f:
        yaml_writer.dump(data, f)

    return defaults


def generate_config_schema(path: str = "config/schema.json") -> dict:
    """Generate JSON schema file and return the schema dict."""
    import json
    import math

    def resolve_refs(obj: dict, defs: dict) -> dict:
        """Resolve nested schema references to show them nicely in the documentation."""
        if isinstance(obj, dict):
            if "$ref" in obj:
                ref_path = obj["$ref"]  # "#/$defs/RunConfig
                ref_name = ref_path.split("/")[-1]
                if ref_name in defs:
                    resolved = resolve_refs(defs[ref_name].copy(), defs)
                    # Keep description from the reference
                    if "description" in obj and "description" not in resolved:
                        resolved["description"] = obj["description"]
                    return resolved
            return {k: resolve_refs(v, defs) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [resolve_refs(item, defs) for item in obj]
        return obj

    def sanitize_for_json(obj):
        """Replace infinity values with None for valid JSON."""
        if isinstance(obj, dict):
            return {k: sanitize_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [sanitize_for_json(v) for v in obj]
        elif isinstance(obj, float) and math.isinf(obj):
            return None
        return obj

    def remove_nested_titles(obj, is_root=True):
        """Remove nested titles (e.g. model class names)."""
        if isinstance(obj, dict):
            result = {}
            for k, v in obj.items():
                if k == "title" and not is_root:
                    continue
                result[k] = remove_nested_titles(v, is_root=False)
            return result
        elif isinstance(obj, list):
            return [remove_nested_titles(item, is_root=False) for item in obj]
        return obj

    def remove_object_type(obj, is_root=True):
        """Remove 'type: object' from nested objects (redundant when properties exist)."""
        if isinstance(obj, dict):
            result = {}
            for k, v in obj.items():
                if (
                    k == "type"
                    and v == "object"
                    and not is_root
                    and "properties" in obj
                ):
                    continue
                result[k] = remove_object_type(v, is_root=False)
            return result
        elif isinstance(obj, list):
            return [remove_object_type(item, is_root=False) for item in obj]
        return obj

    def convert_rst_to_markdown(obj):
        """Convert RST-style links in 'description' to Markdown in 'markdownDescription'."""

        def rst_to_md(text):
            """Convert RST link format `Link Text <URL>`_ to Markdown [Link Text](URL)."""
            # Pattern matches: `Link Text <URL>`_
            pattern = r"`([^<>`]+)\s*<([^>]+)>`_"
            return re.sub(pattern, r"[\1](\2)", text)

        if isinstance(obj, dict):
            result = {}
            for k, v in obj.items():
                if k == "description" and isinstance(v, str) and "`" in v and "<" in v:
                    result[k] = v
                    md_text = rst_to_md(v)
                    if md_text != v:
                        result["markdownDescription"] = md_text
                else:
                    result[k] = convert_rst_to_markdown(v)
            return result
        elif isinstance(obj, list):
            return [convert_rst_to_markdown(item) for item in obj]
        return obj

    schema = ConfigSchema.model_json_schema()
    defs = schema.get("$defs", {})
    schema = resolve_refs(schema, defs)
    schema = sanitize_for_json(schema)
    schema = remove_nested_titles(schema)
    schema = remove_object_type(schema)
    schema = convert_rst_to_markdown(schema)
    with open(path, "w") as f:
        json.dump(schema, f, indent=2)
        f.write("\n")
    return schema


__all__ = [
    "ConfigSchema",
    "validate_config",
    "generate_config_defaults",
    "generate_config_schema",
    "ValidationError",
]
