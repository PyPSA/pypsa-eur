# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Enable configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#enable
"""

from typing import Literal

from pydantic import BaseModel, Field


class EnableConfig(BaseModel):
    """Configuration for `enable` settings."""

    retrieve: bool | Literal["auto"] = Field(
        "auto",
        description="Switch to include (true) or exclude (false) the retrieve_* rules of snakemake into the workflow; 'auto' sets true|false based on availability of an internet connection to prevent issues with snakemake failing due to lack of internet connection.",
    )
    retrieve_databundle: bool = Field(
        True,
        description="Switch to retrieve databundle from zenodo via the rule `retrieve_databundle` or whether to keep a custom databundle located in the corresponding folder.",
    )
    retrieve_cost_data: bool = Field(
        True,
        description="Switch to retrieve technology cost data from [technology-data repository](https://github.com/PyPSA/technology-data).",
    )
    build_cutout: bool = Field(
        False,
        description="Switch to enable the building of cutouts via the rule `build_cutout`.",
    )
    retrieve_cutout: bool = Field(
        True,
        description="Switch to enable the retrieval of cutouts from zenodo with `retrieve_cutout`.",
    )
    drop_leap_day: bool = Field(
        True,
        description="Switch to drop February 29 from all time-dependent data in leap years.",
    )
