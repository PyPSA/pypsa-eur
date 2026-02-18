# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Lines configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#lines
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _DynamicLineRatingConfig(ConfigModel):
    """Configuration for `lines.dynamic_line_rating` settings."""

    activate: bool = Field(
        False,
        description="Whether to take dynamic line rating into account.",
    )
    cutout: str | list[str] = Field(
        "default",
        description="Specifies the weather data cutout file(s) to use.",
    )
    correction_factor: float = Field(
        0.95,
        description="Factor to compensate for overestimation of wind speeds in hourly averaged wind data.",
    )
    max_voltage_difference: float | Literal[False] = Field(
        False,
        description="Maximum voltage angle difference in degrees or 'false' to disable.",
    )
    max_line_rating: float | Literal[False] = Field(
        False,
        description="Maximum line rating relative to nominal capacity without DLR, e.g. 1.3 or 'false' to disable.",
    )


class LinesConfig(BaseModel):
    """Configuration for `lines` settings."""

    types: dict[float, str] = Field(
        default_factory=lambda: {
            63.0: "94-AL1/15-ST1A 20.0",
            66.0: "94-AL1/15-ST1A 20.0",
            90.0: "184-AL1/30-ST1A 110.0",
            110.0: "184-AL1/30-ST1A 110.0",
            132.0: "243-AL1/39-ST1A 110.0",
            150.0: "243-AL1/39-ST1A 110.0",
            220.0: "Al/St 240/40 2-bundle 220.0",
            300.0: "Al/St 240/40 3-bundle 300.0",
            330.0: "Al/St 240/40 3-bundle 300.0",
            380.0: "Al/St 240/40 4-bundle 380.0",
            400.0: "Al/St 240/40 4-bundle 380.0",
            500.0: "Al/St 240/40 4-bundle 380.0",
            750.0: "Al/St 560/50 4-bundle 750.0",
        },
        description="Specifies line types to assume for the different voltage levels of the ENTSO-E grid extraction. Should normally handle voltage levels 220, 300, and 380 kV.",
    )
    s_max_pu: float = Field(
        0.7,
        description="Correction factor for line capacities (`s_nom`) to approximate N-1 security and reserve capacity for reactive power flows.",
    )
    s_nom_max: float = Field(
        float("inf"),
        description="Global upper limit for the maximum capacity of each extendable line (MW).",
    )
    max_extension: float = Field(
        20000,
        description="Upper limit for the extended capacity of each extendable line (MW).",
    )
    length_factor: float = Field(
        1.25,
        description="Correction factor to account for the fact that buses are *not* connected by lines through air-line distance.",
    )
    reconnect_crimea: bool = Field(
        True,
        description="Whether to reconnect Crimea to the Ukrainian grid.",
    )
    under_construction: Literal["zero", "remove", "keep"] = Field(
        "keep",
        description="Specifies how to handle lines which are currently under construction.",
    )
    dynamic_line_rating: _DynamicLineRatingConfig = Field(
        default_factory=_DynamicLineRatingConfig,
        description="Configuration for dynamic line rating.",
    )
