# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Adjustments configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#adjustments
"""

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _AdjustmentConfig(ConfigModel):
    """Configuration for adjustment settings (factor/absolute)"""

    factor: bool | dict[str, dict[str, dict[str, float | dict[int, float]]]] = Field(
        False,
        description="Multiply original value with given factor",
    )
    absolute: bool | dict[str, dict[str, dict[str, float | dict[int, float]]]] = Field(
        False,
        description="Set attribute to absolute value. Can be also a dictionary with planning horizons as keys.",
    )


class AdjustmentsConfig(BaseModel):
    """Configuration for top-level adjustments key."""

    electricity: bool | _AdjustmentConfig = Field(
        False,
        description="Parameter adjustments applied in `prepare_network`.",
    )
    sector: bool | _AdjustmentConfig = Field(
        default_factory=lambda: _AdjustmentConfig(
            factor={
                "Link": {
                    "electricity distribution grid": {
                        "capital_cost": 1.0,
                    }
                }
            },
            absolute=False,
        ),
        description="Parameter adjustments applied in `prepare_sector_network`.",
    )
