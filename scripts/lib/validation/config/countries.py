# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Countries configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#countries
"""

from pydantic import Field, RootModel


class CountriesConfig(RootModel[list[str]]):
    """Configuration for `countries` settings."""

    root: list[str] = Field(
        default=[
            "AL",
            "AT",
            "BA",
            "BE",
            "BG",
            "CH",
            "CZ",
            "DE",
            "DK",
            "EE",
            "ES",
            "FI",
            "FR",
            "GB",
            "GR",
            "HR",
            "HU",
            "IE",
            "IT",
            "LT",
            "LU",
            "LV",
            "ME",
            "MK",
            "NL",
            "NO",
            "PL",
            "PT",
            "RO",
            "RS",
            "SE",
            "SI",
            "SK",
            "XK",
        ],
        description="European countries defined by their `Two-letter country codes (ISO 3166-1) <https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2>`_ which should be included in the energy system model.",
    )
