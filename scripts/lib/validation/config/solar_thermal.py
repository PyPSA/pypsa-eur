# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Solar thermal configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#solar-thermal
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _OrientationConfig(ConfigModel):
    """Configuration for `solar_thermal.orientation` settings."""

    slope: float = Field(
        45.0,
        description="The angle between the ground and the panels.",
    )
    azimuth: float = Field(
        180.0,
        description="The angle between the North and the sun with panels on the local horizon.",
    )


class SolarThermalConfig(BaseModel):
    """Configuration for `solar_thermal` settings."""

    clearsky_model: Literal["simple", "enhanced"] = Field(
        "simple",
        description="Type of clearsky model for diffuse irradiation.",
    )
    orientation: _OrientationConfig = Field(
        default_factory=_OrientationConfig,
        description="Panel orientation with slope and azimuth.",
    )
    cutout: str = Field(
        "default",
        description="Name of the cutout to use for solar thermal calculations.",
    )
