# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Conventional generators configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration/#conventional_cf
"""

from pydantic import BaseModel, ConfigDict, Field

from scripts.lib.validation.config._base import ConfigModel


class _EfficiencyParameters(BaseModel):
    """Parameters of the linear efficiency heuristic for one carrier."""

    efficiency: float = Field(description="Efficiency of a plant built in `year`.")
    slope: float = Field(description="Efficiency gain per later build year.")
    year: int = Field(description="Reference build year of `efficiency`.")
    max: float = Field(description="Upper bound of the estimated efficiency.")


class _EstimateEfficienciesConfig(BaseModel):
    """Configuration for `conventional.estimate_efficiencies` settings."""

    enable: bool = Field(
        False,
        description="Estimate missing plant-level efficiencies from a carrier- and age-dependent linear heuristic. The build year is the year of the last retrofit or, if unknown, of commissioning. Plants without build year or carriers without parameters fall back to the technology cost data.",
    )
    reference_year: int = Field(
        2025,
        description="Year in which the age of power plants is evaluated for degradation.",
    )
    degradation_start: int = Field(
        10,
        description="Age in years after which efficiency degradation starts.",
    )
    degradation_rate: float = Field(
        0.001,
        description="Relative efficiency loss per year of age after `degradation_start`.",
    )
    parameters: dict[str, _EfficiencyParameters] = Field(
        default_factory=lambda: {
            carrier: _EfficiencyParameters(efficiency=e, slope=s, year=y, max=m)
            for carrier, (e, s, y, m) in {
                "lignite": (0.25, 0.003, 1960, 0.42),
                "coal": (0.28, 0.003, 1960, 0.44),
                "CCGT": (0.40, 0.004, 1980, 0.60),
                "OCGT": (0.28, 0.003, 1970, 0.41),
                "oil": (0.28, 0.002, 1960, 0.38),
                "nuclear": (0.33, 0.0, 1960, 0.33),
            }.items()
        },
        description="Per carrier: efficiency `efficiency` of a plant built in `year`, rising by `slope` per later build year up to `max`.",
    )


class ConventionalConfig(ConfigModel):
    """Configuration for `conventional` settings."""

    model_config = ConfigDict(extra="allow")

    estimate_efficiencies: _EstimateEfficienciesConfig = Field(
        default_factory=_EstimateEfficienciesConfig,
        description="Estimation of missing plant-level efficiencies.",
    )
    unit_commitment: bool = Field(
        False,
        description="Allow the overwrite of ramp_limit_up, ramp_limit_start_up, ramp_limit_shut_down, p_min_pu, min_up_time, min_down_time, and start_up_cost of conventional generators. Refer to the CSV file 'unit_commitment.csv'.",
    )
    dynamic_fuel_price: bool = Field(
        False,
        description="Consider the monthly fluctuating fuel prices for each conventional generator. Refer to the CSV file 'data/validation/monthly_fuel_price.csv'.",
    )
    fuel_price_rolling_window: int = Field(
        6,
        description="Monthly rolling mean window for fossil fuel prices smoothing.",
        ge=1,
    )
    nuclear: dict[str, str | float] = Field(
        default_factory=lambda: {"p_max_pu": "data/nuclear_p_max_pu.csv"},
        description="For any carrier/technology overwrite attributes as listed below.",
    )
