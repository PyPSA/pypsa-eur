# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Transmission candidate configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#transmission
"""

from typing import Literal

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


class _TransmissionCarrierConfigGeneral(BaseModel):
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


class _TransmissionCarrierConfigElectricity(BaseModel):
    """Configuration for electricity transmission grid."""

    enable: bool = Field(
        True,
        description="Switch for enabling/disabling the electricity transmission grid.",
    )
    base_network: Literal["entsoegridkit", "osm", "tyndp"] = Field(
        "osm",
        description="Specify the underlying base network, i.e. GridKit (based on ENTSO-E web map extract), OpenStreetMap (OSM), or TYNDP.",
    )
    transmission_limit: str = Field(
        "vopt",
        description="Limit on transmission expansion. The first part can be `v` (for setting a limit on line volume) or `c` (for setting a limit on line cost). The second part can be `opt` or a float bigger than one (e.g. 1.25). If `opt` is chosen line expansion is optimised according to its capital cost (where the choice `v` only considers overhead costs for HVDC transmission lines, while `c` uses more accurate costs distinguishing between overhead and underwater sections and including inverter pairs). The setting `v1.25` will limit the total volume of line expansion to 25% of currently installed capacities weighted by individual line lengths. The setting `c1.25` will allow to build a transmission network that costs no more than 25 % more than the current system.",
    )


class TransmissionConfig(BaseModel):
    """Configuration for `transmission` settings."""

    electricity: _TransmissionCarrierConfigElectricity = Field(
        default_factory=_TransmissionCarrierConfigElectricity,
        description="Configuration for electricity transmission grid.",
    )
    hydrogen: _TransmissionCarrierConfigGeneral = Field(
        default_factory=_TransmissionCarrierConfigGeneral,
        description="Configuration for hydrogen transmission candidates.",
    )
    carbon_dioxide: _TransmissionCarrierConfigGeneral = Field(
        default_factory=_TransmissionCarrierConfigGeneral,
        description="Configuration for carbon dioxide transmission candidates.",
    )
