# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Scenario configuration block.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#scenario
Wildcard docs in https://pypsa-eur.readthedocs.io/en/latest/wildcards.html
"""

from pydantic import BaseModel, Field


class ScenarioConfig(BaseModel):
    """Configuration for top level `scenario` settings."""

    clusters: list[int] = Field(
        default_factory=lambda: [50],
        description="List of ``{clusters}`` wildcards to run.",
    )
    opts: list[str] = Field(
        default_factory=lambda: [""],
        description="List of ``{opts}`` wildcards to run.",
    )
    sector_opts: list[str] = Field(
        default_factory=lambda: [""],
        description="List of ``{sector_opts}`` wildcards to run.",
    )
    planning_horizons: list[int] = Field(
        default_factory=lambda: [2050],
        description="List of ``{planning_horizon}`` wildcards to run.",
    )
