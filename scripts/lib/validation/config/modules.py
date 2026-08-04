# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Modules configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#modules
"""

from pathlib import Path

from pydantic import Field, FilePath

from scripts.lib.validation.config._base import ConfigModel


class _ModuleConfig(ConfigModel):
    """Configuration for module settings"""

    config_path: FilePath = Field(
        description="Path to module default configuration file. Can be relative to the project directory or absolute.",
    )
    version: str = Field(
        description="Module version to use.",
    )


class _GeoBoundariesModuleConfig(_ModuleConfig):
    """Configuration for module settings"""

    scenario: str = Field(
        default="default",
        description="Scenario to use for the geo_boundaries module.",
    )


class ModulesConfig(ConfigModel):
    """Configuration for modules."""

    geo_boundaries: _GeoBoundariesModuleConfig = Field(
        default=_GeoBoundariesModuleConfig(
            config_path=Path("config/modules/geo_boundaries.yaml"), version="v1.0.1"
        ),
        description="Configuration for the geo_boundaries module.",
    )
