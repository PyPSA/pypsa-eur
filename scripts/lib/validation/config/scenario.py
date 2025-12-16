# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Scenario configuration block.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#scenario
Wildcard docs in https://pypsa-eur.readthedocs.io/en/latest/wildcards.html
"""

from typing import Literal

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class ScenarioConfig(ConfigModel):
    """Configuration for top level `scenario` settings."""

    clusters: list[int | Literal["adm", "all"]] = Field(
        default_factory=lambda: [50],
        description="List of ``{clusters}`` wildcards to run. Use 'adm' for administrative clustering mode, 'all' for all nodes.",
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
