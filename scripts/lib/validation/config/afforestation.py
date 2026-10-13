# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Afforestation configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#afforestation
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class AfforestationConfig(ConfigModel):
    """Configuration for `afforestation` settings (used with `sector: afforestation: true`)."""

    use_discount_rate: bool = Field(
        False,
        description="Annualise the investment with the discount rate of the cost data; if false, without discounting.",
    )
    crcf_efficiency: float = Field(
        0.8,
        description="Share of the gross sequestration credited as net removal under the EU Carbon Removals and Carbon Farming regulation (EU 2024/3012): uncertainty deduction (at least 10%) and buffer-pool contribution (8-17% risk rate), (1-0.10)*(1-0.11) = 0.80, range 0.72-0.83.",
    )
