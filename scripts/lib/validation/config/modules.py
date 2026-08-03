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

    default_config: FilePath = Field(
        description="Path to module default configuration file. Can be relative to the project directory or absolute.",
    )
    version: str = Field(
        description="Module version to use.",
    )
    config: dict = Field(
        default_factory=dict,
        description="Module configuration overrides. Overrides the default configuration.",
    )


class ModulesConfig(ConfigModel):
    """Configuration for modules."""

    geo_boundaries: _ModuleConfig = Field(
        default=_ModuleConfig(
            default_config=Path("config/modules/geo_boundaries.yaml"),
            version="v1.0.1",
        ),
        description="Configuration for the geo_boundaries module.",
    )
