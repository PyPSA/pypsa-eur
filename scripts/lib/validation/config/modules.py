# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Modules configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration/#modules_cf
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class _GeoBoundariesModuleConfig(ConfigModel):
    """Configuration for the geo_boundaries module."""

    version: str = Field(
        default="v1.0.1",
        description="Release tag of the geo_boundaries module. Applies to all scenarios of a run.",
    )
    nuts_year: int = Field(
        default=2021,
        description="NUTS release year used for countries with NUTS3 regions.",
    )
    nuts_resolution: str = Field(
        default="01M",
        pattern=r"^[0-9]{2}M$",
        description="Resolution of the NUTS shapes, e.g. '01M' for 1:1 million.",
    )
    adm1_countries: list[str] = Field(
        default=["BA", "MD", "UA", "XK"],
        description="Countries without NUTS3 regions. Their ADM1 regions are taken from geoBoundaries (gbOpen release).",
    )


class ModulesConfig(ConfigModel):
    """Configuration for external Modelblocks modules."""

    geo_boundaries: _GeoBoundariesModuleConfig = Field(
        default_factory=_GeoBoundariesModuleConfig,
        description="Configuration for the geo_boundaries module.",
    )
