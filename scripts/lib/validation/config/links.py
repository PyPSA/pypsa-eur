# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Links configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#links
"""

from typing import Literal

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class LinksConfig(ConfigModel):
    """Configuration for `links` settings."""

    p_max_pu: float = Field(
        1.0,
        description="Correction factor for link capacities `p_nom`.",
    )
    p_min_pu: float = Field(
        -1.0,
        description="Correction factor for link capacities `p_nom`.",
    )
    p_nom_max: float = Field(
        float("inf"),
        description="Global upper limit for the maximum capacity of each extendable DC link (MW).",
    )
    max_extension: float = Field(
        30000,
        description="Upper limit for the extended capacity of each extendable DC link (MW).",
    )
    length_factor: float = Field(
        1.25,
        description="Correction factor to account for the fact that buses are *not* connected by links through air-line distance.",
    )
    under_construction: Literal["zero", "remove", "keep"] = Field(
        "keep",
        description="Specifies how to handle lines which are currently under construction.",
    )
