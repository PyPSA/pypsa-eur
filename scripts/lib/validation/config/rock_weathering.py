# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Enhanced rock weathering configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#rock_weathering
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class RockWeatheringConfig(ConfigModel):
    """Configuration for `rock_weathering` settings."""

    co2_removal_per_sqkm: float = Field(
        309,
        description="Tonnes of CO2 sequestered per square kilometre of eligible land per year (temperate-only rate).",
    )
    max_land_usage: float = Field(
        0.2,
        description="Maximum fraction of suitable area used for rock weathering.",
    )
