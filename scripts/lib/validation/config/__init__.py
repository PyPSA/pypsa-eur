# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Config validation for PyPSA-EUR.
"""

from pydantic import Field, ValidationError

from scripts.lib.validation.config.run import RunConfig
from scripts.lib.validation.config.top_level import TopLevelConfig


class ConfigSchema(TopLevelConfig):
    """
    Combined configuration schema for PyPSA-EUR.
    """

    run: RunConfig = Field(
        default_factory=RunConfig,
        description="Run configuration for PyPSA-EUR workflow execution.",
    )


def validate_config(config: dict) -> ConfigSchema:
    """Validate config dict against schema."""
    return ConfigSchema(**config)


def export_defaults() -> dict:
    """Export default config values as dict."""
    return ConfigSchema().model_dump()


def export_json_schema() -> dict:
    """Export JSON schema for ConfigSchema."""
    return ConfigSchema.model_json_schema()


__all__ = [
    "ConfigSchema",
    "validate_config",
    "export_json_schema",
    "export_defaults",
    "ValidationError",
]
