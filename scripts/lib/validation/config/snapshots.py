# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Snapshots configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#snapshots
"""

from typing import Literal

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class SnapshotsConfig(ConfigModel):
    """Configuration for `snapshots` settings."""

    start: str | list[str] = Field(
        "2013-01-01",
        description="Left bound of date range.",
    )
    end: str | list[str] = Field(
        "2014-01-01",
        description="Right bound of date range.",
    )
    inclusive: Literal["left", "right", "both"] | None = Field(
        "left",
        description="Make the time interval closed to the `left`, `right`, or both sides `both` or neither side `None`.",
    )
