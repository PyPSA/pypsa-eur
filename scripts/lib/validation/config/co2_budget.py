# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
CO2 budget configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#co2-budget
"""

from pydantic import Field, RootModel


class Co2BudgetConfig(RootModel[dict[int, float]]):
    """Configuration for `co2_budget` settings."""

    root: dict[int, float] = Field(
        default_factory=lambda: {
            # CO2 budget values as fraction of 1990 emissions per planning horizon year.
            2020: 0.720,  # Average emissions of 2019-2021 relative to 1990, CO2 excl LULUCF
            2025: 0.648,  # WAM projection from EEA Member States' GHG projections 2023
            2030: 0.450,  # 55% reduction by 2030 (Fit for 55)
            2035: 0.250,  # Interpolated
            2040: 0.100,  # 90% reduction target by 2040
            2045: 0.050,  # Interpolated
            2050: 0.000,  # Climate-neutral by 2050
            # Sources:
            # - EEA Annual European Union GHG inventory 1990-2021 (https://unfccc.int/documents/627830)
            # - EEA Member States' GHG projections 2023 (https://www.eea.europa.eu/en/datahub/datahubitem-view/4b8d94a4-aed7-4e67-a54c-0623a50f48e8)
        },
        description="CO2 budget as a fraction of 1990 emissions. Overwritten if `Co2Lx` or `cb` are set in `{sector_opts}` wildcard.",
    )
