# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Load configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#load
"""

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _FillGapsConfig(ConfigModel):
    """Configuration for `load.fill_gaps` settings."""

    enable: bool = Field(
        True,
        description="Whether to fill gaps using interpolation for small gaps and time shift for large gaps.",
    )
    interpolate_limit: int = Field(
        3,
        description="Maximum gap size (consecutive nans) which interpolated linearly.",
    )
    time_shift_for_large_gaps: str = Field(
        "1w",
        description="Periods which are used for copying time-slices in order to fill large gaps of nans. Have to be valid `pandas` period strings.",
    )


class _DistributionKeyConfig(BaseModel):
    """Configuration for `load.distribution_key` settings."""

    gdp: float = Field(
        0.6,
        description="Weighting factor for the GDP data in the distribution key.",
    )
    population: float = Field(
        0.4,
        description="Weighting factor for the population data in the distribution key.",
    )


class LoadConfig(BaseModel):
    """Configuration for `load` settings."""

    fill_gaps: _FillGapsConfig = Field(
        default_factory=_FillGapsConfig,
        description="Gaps filling strategy used.",
    )
    manual_adjustments: bool = Field(
        True,
        description="Whether to adjust the load data manually according to the function in `manual_adjustment`.",
    )
    scaling_factor: float = Field(
        1.0,
        description="Global correction factor for the load time series.",
    )
    fixed_year: int | bool = Field(
        False,
        description="To specify a fixed year for the load time series that deviates from the snapshots' year.",
    )
    supplement_synthetic: bool = Field(
        True,
        description="Whether to supplement missing data for selected time period should be supplemented by synthetic data from `Zenodo <https://zenodo.org/records/10820928>`_.",
    )
    distribution_key: _DistributionKeyConfig = Field(
        default_factory=_DistributionKeyConfig,
        description="Distribution key for spatially disaggregating the per-country electricity demand data.",
    )
