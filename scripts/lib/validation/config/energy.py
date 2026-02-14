# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Energy configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#energy
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class EnergyConfig(ConfigModel):
    """Configuration for `energy` settings."""

    energy_totals_year: int = Field(
        2023,
        description="The year for the sector energy use. The year must be available in the Eurostat report.",
    )
    base_emissions_year: int = Field(
        1990,
        description="The base year for the sector emissions. See `European Environment Agency (EEA) <https://www.eea.europa.eu/data-and-maps/data/national-emissions-reported-to-the-unfccc-and-to-the-eu-greenhouse-gas-monitoring-mechanism-16>`_.",
    )
    emissions: str = Field(
        "CO2",
        description="Specify which sectoral emissions are taken into account. Data derived from EEA. Currently only CO2 is implemented.",
    )
