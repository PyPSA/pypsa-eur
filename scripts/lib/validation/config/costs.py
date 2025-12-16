# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Costs configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#costs
"""

from pydantic import BaseModel, ConfigDict, Field

from scripts.lib.validation.config._base import ConfigModel


class _EmissionPricesConfig(ConfigModel):
    """Configuration for `costs.emission_prices` settings."""

    enable: bool = Field(
        False,
        description="Add cost for a carbon-dioxide price configured in `costs: emission_prices: co2` to `marginal_cost` of generators. Config setting can also be enabled with the keyword `Ep` in the `{opts}` wildcard for electricity-only runs.",
    )
    co2: float | dict[str, float] = Field(
        0.0,
        description="Exogenous price of carbon-dioxide. In electricity-only runs it is added to the marginal costs of fossil-fuelled generators according to their carbon intensity, while for sector networks it applies to emissions ending up in CO2 atmosphere.",
    )
    co2_monthly_prices: bool = Field(
        False,
        description="Add monthly cost for a carbon-dioxide price based on historical values built by the rule `build_monthly_prices`.",
    )


class _FillValuesConfig(BaseModel):
    """Configuration for `costs.fill_values` settings."""

    FOM: float = Field(0, description="Default fixed operation and maintenance cost.")
    VOM: float = Field(
        0, description="Default variable operation and maintenance cost."
    )
    efficiency: float = Field(1, description="Default efficiency.")
    fuel: float = Field(0, description="Default fuel cost.")
    investment: float = Field(0, description="Default investment cost.")
    lifetime: int = Field(25, description="Default lifetime in years.")
    CO2_intensity: float = Field(
        0, alias="CO2 intensity", description="Default CO2 intensity."
    )
    discount_rate: float = Field(
        0.07, alias="discount rate", description="Default discount rate."
    )
    standing_losses: float = Field(
        0, alias="standing losses", description="Default standing losses."
    )

    model_config = ConfigDict(populate_by_name=True)


class CostsConfig(BaseModel):
    """Configuration for `costs` settings."""

    year: int = Field(
        2050,
        description="Year for which to retrieve cost assumptions of `data/costs/primary/<version>/costs_<year>.csv`.",
    )
    social_discountrate: float = Field(
        0.02,
        description="Social discount rate to compare costs in different investment periods. 0.02 corresponds to a social discount rate of 2%.",
    )
    fill_values: _FillValuesConfig = Field(
        default_factory=_FillValuesConfig,
        description="Default values if not specified for a technology in `resources/costs.csv`.",
    )
    custom_cost_fn: str | None = Field(
        "data/custom_costs.csv",
        description="Path to the custom costs file. None if it should not be used. Default `data/custom_costs.csv` contains minor adjustments for stabilising the optimisation results.",
    )
    overwrites: dict[str, dict[str, float]] = Field(
        default_factory=dict,
        description="For the given parameters and technologies, assumptions about their parameter are overwritten the corresponding value of the technology.",
    )
    capital_cost: dict[str, float] = Field(
        default_factory=dict,
        description="For the given technologies, assumptions about their capital investment costs are set to the corresponding value. Optional; overwrites cost assumptions from `resources/costs.csv`.",
    )
    marginal_cost: dict[str, float] = Field(
        default_factory=dict,
        description="For the given technologies, assumptions about their marginal operating costs are set to the corresponding value. Optional; overwrites cost assumptions from `resources/costs.csv`.",
    )
    emission_prices: _EmissionPricesConfig = Field(
        default_factory=_EmissionPricesConfig,
        description="Specify exogenous prices for emission types listed in `network.carriers` to marginal costs.",
    )
