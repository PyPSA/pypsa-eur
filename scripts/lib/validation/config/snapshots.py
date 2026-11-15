# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Snapshots configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration/#snapshots_cf
"""

import re
from typing import Any

from pydantic import Field, model_validator

from scripts.lib.validation.config._base import ConfigModel

MIGRATION_HINT = (
    "`snapshots: end` is now the last snapshot and is included, as in "
    "`pandas.date_range`, and `snapshots: inclusive` was removed. "
    'Replace e.g. `end: "2014-01-01"` with `end: "2013-12-31 23:00"`. '
    "See the release notes for a migration guide."
)


class SnapshotsConfig(ConfigModel):
    """Configuration for `snapshots` settings."""

    start: str | list[str] = Field(
        "2013-01-01 00:00",
        description="First snapshot (included).",
    )
    end: str | list[str] = Field(
        "2013-12-31 23:00",
        description="Last snapshot (included), with time, e.g. `2013-12-31 23:00`.",
    )

    @model_validator(mode="before")
    @classmethod
    def reject_legacy_settings(cls, data: Any) -> Any:
        """Reject `inclusive` and a date-only `end`, which was exclusive before."""
        if not isinstance(data, dict):
            return data
        if "inclusive" in data:
            raise ValueError(MIGRATION_HINT)
        ends = data.get("end", [])
        for end in ends if isinstance(ends, list) else [ends]:
            if re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(end)):
                raise ValueError(
                    f"`snapshots: end: {end}` needs a time, e.g. "
                    f'`end: "{end} 00:00"`. {MIGRATION_HINT}'
                )
        return data
