# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
CO2 budget configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#co2-budget
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class Co2BudgetConfig(ConfigModel):
    """Configuration for `co2_budget` settings."""

    emissions_scope: str = Field(
        "All greenhouse gases - (CO2 equivalent)",
        description="Emissions scope for CO2 budget calculations.",
    )
    budget_type: str = Field(
        "fraction",
        alias="values",
        description="Interpretation of budget values: 'fraction' or 'absolute'.",
    )
    upper: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.72,
            2025: 0.648,
            2030: 0.45,
            2035: 0.25,
            2040: 0.1,
            2045: 0.05,
            2050: 0.0,
        },
        description="Upper CO2 budget as fraction of base year emissions per planning horizon.",
    )
    lower: dict[int, float] | None = Field(
        default=None,
        description="Lower CO2 budget limit to prevent unrealistic/too-fast decarbonization.",
    )
