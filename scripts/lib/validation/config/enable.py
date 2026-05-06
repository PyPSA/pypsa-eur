# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Enable configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#enable
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class EnableConfig(ConfigModel):
    """Configuration for `enable` settings."""

    drop_leap_day: bool = Field(
        True,
        description="Switch to drop February 29 from all time-dependent data in leap years.",
    )
