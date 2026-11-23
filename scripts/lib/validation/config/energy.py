# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Energy configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration/#energy_cf
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class EnergyConfig(ConfigModel):
    """Configuration for `energy` settings."""

    energy_totals_year: int = Field(
        2023,
        description="The year for the sector energy use. The year must be available in the Eurostat report.",
    )
