# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Transmission candidate configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#transmission
"""

from pydantic import BaseModel, Field


class _GabrielFilterConfig(BaseModel):
    """Configuration for `transmission.<carrier>.gabriel_filter` settings."""

    enable: bool = Field(
        True,
        description="Whether to filter Delaunay edges to Gabriel edges before min-degree backfilling.",
    )
    min_degree: int = Field(
        1,
        ge=0,
        description="Minimum node degree target applied after Gabriel filtering.",
    )


class _TransmissionCarrierConfig(BaseModel):
    """Configuration for a single transmission carrier."""

    enable: bool = Field(
        True,
        description="Enable transmission candidate generation for this carrier.",
    )
    gabriel_filter: _GabrielFilterConfig = Field(
        default_factory=_GabrielFilterConfig,
        description="Gabriel filter configuration.",
    )
    length_factor: float = Field(
        1.25,
        gt=0,
        description="Multiplier applied to geometric corridor lengths.",
    )
    cost_factor: float = Field(
        1,
        gt=0,
        description="Multiplier applied to the capital cost of transmission infrastructure.",
    )


class TransmissionConfig(BaseModel):
    """Configuration for `transmission` settings."""

    hydrogen: _TransmissionCarrierConfig = Field(
        default_factory=_TransmissionCarrierConfig,
        description="Configuration for hydrogen transmission candidates.",
    )
    carbon_dioxide: _TransmissionCarrierConfig = Field(
        default_factory=_TransmissionCarrierConfig,
        description="Configuration for carbon dioxide transmission candidates.",
    )
