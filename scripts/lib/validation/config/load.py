# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Load configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration/#load_cf
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _FillGapsConfig(ConfigModel):
    """Configuration for `load.fill_gaps` settings."""

    enable: bool = Field(
        True,
        description="Whether to fill gaps using interpolation for small gaps and time shift for large gaps.",
    )
    interpolate_limit: int = Field(
        6,
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


class _DeviceDRConfig(ConfigModel):
    """One device's demand-response settings."""

    shift_time: int = Field(
        ...,
        ge=1,
        description="Hours within which shifted load for this device must be recovered (S).",
    )
    direction: Literal["delay", "advance", "both"] = Field(
        ...,
        description="Which way this device's load can be shifted.",
    )


def _default_dr_devices() -> dict[str, _DeviceDRConfig]:
    """
    Default demand-response devices: shift times and directions from
    O'Reilly et al. (2022)'s own "energy shift constraints" table for the
    residential DR dataset (Zenodo 10.5281/zenodo.7871106), not assumed.
    """
    return {
        "WM": _DeviceDRConfig(shift_time=6, direction="both"),  # washing machine
        "TD": _DeviceDRConfig(shift_time=6, direction="both"),  # tumble drier
        "DW": _DeviceDRConfig(shift_time=6, direction="both"),  # dish washer
        "EV": _DeviceDRConfig(shift_time=4, direction="delay"),  # battery electric vehicle
        "RF": _DeviceDRConfig(shift_time=2, direction="delay"),  # refrigerator
        "FR": _DeviceDRConfig(shift_time=2, direction="delay"),  # freezer
        "AC": _DeviceDRConfig(shift_time=2, direction="delay"),  # air conditioning
        "CP": _DeviceDRConfig(shift_time=2, direction="delay"),  # heat circulation pump
        "HP": _DeviceDRConfig(shift_time=2, direction="advance"),  # air-to-air heat pump
        "SH": _DeviceDRConfig(shift_time=12, direction="advance"),  # storage heater
        "WH": _DeviceDRConfig(shift_time=12, direction="advance"),  # water heater
    }


class _DemandResponseConfig(ConfigModel):
    """Configuration for `load.demand_response` settings."""

    enable: bool = Field(
        False, description="Whether to add residential demand response."
    )
    data_dir: str = Field(
        "data/demand_response/achievable",
        description="Directory with prepare_dr_data.py's per-device, per-year CSV output.",
    )
    devices: dict[str, _DeviceDRConfig] = Field(
        default_factory=_default_dr_devices,
        description="Devices to model, keyed by O'Reilly device code (WM, TD, DW, EV, RF, FR, AC, CP, HP, SH, WH). "
        "Default shift times and directions are from O'Reilly et al. (2022)'s own "
        '"energy shift constraints" table for the residential DR dataset.',
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
    substation_only: bool = Field(
        True,
        description="Whether to only consider substations for the spatial disaggregation of the per-country electricity demand data.",
    )
    distribution_key: _DistributionKeyConfig = Field(
        default_factory=_DistributionKeyConfig,
        description="Distribution key for spatially disaggregating the per-country electricity demand data.",
    )
    demand_response: _DemandResponseConfig = Field(
        default_factory=_DemandResponseConfig,
        description="Residential demand response settings.",
    )
