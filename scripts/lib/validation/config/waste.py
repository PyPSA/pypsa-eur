# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Waste configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration/#waste_cf
"""

from pydantic import BaseModel, Field


class WasteConfig(BaseModel):
    """Configuration for `waste` settings."""

    transport: bool = Field(
        False,
        description="Switch to enable the transport of waste between bus regions. Only works if `sector: waste_spatial` is activated.",
    )
    waste_to_energy: bool = Field(
        False,
        description="Switch to enable expansion of waste to energy CHPs for conversion of plastics.",
    )
    waste_to_energy_cc: bool = Field(
        False,
        description="Switch to enable expansion of waste to energy CHPs for conversion of plastics with carbon capture.",
    )
