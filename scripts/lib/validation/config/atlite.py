# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Atlite configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#atlite
"""

from pydantic import BaseModel, Field, field_validator

from scripts.lib.validation.config._base import ConfigModel


class _CutoutConfig(ConfigModel):
    """Configuration for a single cutout in `atlite.cutouts`."""

    module: str | list[str] = Field(
        ...,
        description="Source of the reanalysis weather dataset (e.g. `ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ or `SARAH-3 <https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V002>`_).",
    )
    x: list[float] | None = Field(
        None,
        description="Range of longitudes [°] to download weather data for. Float interval within [-180, 180]. If not defined, it defaults to the spatial bounds of all bus shapes.",
    )
    y: list[float] | None = Field(
        None,
        description="Range of latitudes [°] to download weather data for. Float interval within [-90, 90]. If not defined, it defaults to the spatial bounds of all bus shapes.",
    )
    dx: float | None = Field(
        None,
        gt=0.25,
        description="Grid resolution [°] for longitude. Must be larger than 0.25°.",
    )
    dy: float | None = Field(
        None,
        gt=0.25,
        description="Grid resolution [°] for latitude. Must be larger than 0.25°.",
    )
    time: list[str] | None = Field(
        None,
        description="Time span to download weather data for. If not defined, it defaults to the time interval spanned by the snapshots.",
    )
    features: str | list[str] | None = Field(
        None,
        description="When freshly building a cutout, retrieve data only for those features. If not defined, it defaults to all available features.",
    )
    sarah_dir: str | None = Field(
        None,
        description="Path to the location where SARAH-2 or SARAH-3 data is stored; SARAH data requires a manual separate download, see the https://atlite.readthedocs.io for details. Required for building cutouts with SARAH, not required for ERA5 cutouts.",
    )

    @field_validator("x")
    @classmethod
    def validate_longitude(cls, v):
        if v is not None:
            if len(v) != 2:
                raise ValueError("x must be a list of two floats [min, max]")
            if not all(-180 <= val <= 180 for val in v):
                raise ValueError("Longitude values must be within [-180, 180]")
            if v[0] >= v[1]:
                raise ValueError("x[0] must be less than x[1]")
        return v

    @field_validator("y")
    @classmethod
    def validate_latitude(cls, v):
        if v is not None:
            if len(v) != 2:
                raise ValueError("y must be a list of two floats [min, max]")
            if not all(-90 <= val <= 90 for val in v):
                raise ValueError("Latitude values must be within [-90, 90]")
            if v[0] >= v[1]:
                raise ValueError("y[0] must be less than y[1]")
        return v


class AtliteConfig(BaseModel):
    """Configuration for `atlite` settings."""

    cutout_directory: str = Field(
        "cutouts",
        description="Directory to store cutouts.",
    )
    default_cutout: str | list[str] = Field(
        "europe-2013-sarah3-era5",
        description="Defines a default cutout. Can refer to a single cutout or a list of cutouts.",
    )
    nprocesses: int = Field(
        16,
        description="Number of parallel processes in cutout preparation.",
    )
    show_progress: bool = Field(
        False,
        description="Whether progressbar for atlite conversion processes should be shown. False saves time.",
    )
    cutouts: dict[str, _CutoutConfig] = Field(
        default_factory=lambda: {
            "europe-2013-sarah3-era5": _CutoutConfig(
                module=["sarah", "era5"],
                x=[-12.0, 42.0],
                y=[33.0, 72.0],
                dx=0.3,
                dy=0.3,
                time=["2013", "2013"],
            )
        },
        description="Named cutout configurations.",
    )
