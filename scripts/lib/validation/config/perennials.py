# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Perennialisation configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#perennials
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class PerennialsConfig(ConfigModel):
    """Configuration for `perennials` settings."""

    sequestration_co2: float = Field(
        2,
        description="Tonnes of CO2 equivalent sequestered per hectare per year when 1st-generation biofuel cropland is converted to perennial grasses.",
    )
