# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
PyPSA-Eur component configuration.

Regulates what components with which carriers are kept from PyPSA-Eur.
Some technologies are removed because they are implemented differently
(e.g. battery or H2 storage) or have different year-dependent costs.
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class PypsaEurConfig(ConfigModel):
    """Configuration for `pypsa_eur` settings."""

    Bus: list[str] = Field(
        default_factory=lambda: ["AC"],
        description="Bus carriers to keep from PyPSA-Eur.",
    )
    Link: list[str] = Field(
        default_factory=lambda: ["DC"],
        description="Link carriers to keep from PyPSA-Eur.",
    )
    Generator: list[str] = Field(
        default_factory=lambda: [
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "offwind-float",
            "solar-hsat",
            "solar",
            "ror",
            "nuclear",
        ],
        description="Generator carriers to keep from PyPSA-Eur.",
    )
    StorageUnit: list[str] = Field(
        default_factory=lambda: ["PHS", "hydro"],
        description="StorageUnit carriers to keep from PyPSA-Eur.",
    )
    Store: list[str] = Field(
        default_factory=list,
        description="Store carriers to keep from PyPSA-Eur.",
    )
