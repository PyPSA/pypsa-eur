# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Perennialisation configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#perennials
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class PerennialsConfig(ConfigModel):
    """Configuration for `perennials` settings."""

    co2_removal: float = Field(
        2.2,
        description="CO2 removed in the soil (tCO2/ha/yr) when 1st-generation biofuel cropland is converted to grass-clover: soil carbon gain of 0.6 tC/ha/yr (DCA report 193, Table 3.4).",
    )
    co2_removal_efficiency: float = Field(
        0.7,
        description="Credited share of `co2_removal`, accounting for uncertainty, saturation and reversal of soil carbon (carbon farming certification methodologies under Regulation (EU) 2024/3012, July 2026: uncertainty deduction of at least 10% and liability for reversals).",
    )
    emissions_avoided: float = Field(
        1.1,
        description="Field N2O avoided (tCO2e/ha/yr) when unfertilised grass-clover replaces a 1st-generation biofuel crop: maize minus grass-clover (DCA report 193, Tables 3.3 and 3.6).",
    )
    emissions_avoided_efficiency: float = Field(
        0.9,
        description="Credited share of `emissions_avoided`: minimum uncertainty deduction of 10% for soil emission reductions (carbon farming certification methodologies under Regulation (EU) 2024/3012, July 2026).",
    )
    biofuels_1G_land_year: int | None = Field(
        2025,
        description="Year of the JRC ENSPRESO 1st-generation biofuel crop potentials (interpolated between the ENSPRESO years) that sets the cropland available for perennialisation in every planning horizon, phased in with `biomass: share_sustainable_potential_available`. If null, the cropland follows the ENSPRESO potentials of each planning horizon.",
    )
