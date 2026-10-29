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

    sequestration_co2_credited: float = Field(
        2,
        description="CO2 removal credited as carbon dioxide removal (tCO2/ha/yr) when 1st-generation biofuel cropland is converted to grass-clover. It includes the credited share of the soil carbon gain, e.g. 2.2 tCO2/ha/yr for grass-clover (0.6 tC/ha/yr, DCA report 193, Table 3.4, at any N level) minus the uncertainty deduction and the reversals at ploughing when the sward is re-sown every 4-8 years (carbon farming certification methodologies under Regulation (EU) 2024/3012). Grass-clover is fertilised like the 1st-generation crops it replaces, so the field N2O and the ammonia demand do not change.",
    )
    land_multiplier: float = Field(
        1.0,
        ge=0,
        description="User-driven multiplier of the land available for perennials, e.g. to project it or to reach a given potential. 1: the 1st-generation biofuel cropland of ENSPRESO 2020, scaled to the fuel energy of the Eurostat unsustainable bioliquids. Above 1, perennials may also take cropland that grows food or feed (indirect land use not modelled); the crop mix and the distribution over regions do not change.",
    )
