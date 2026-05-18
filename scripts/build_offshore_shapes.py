# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Creates GIS shape files of offshore exclusive economic zones (EEZ).
"""

import logging

import country_converter as coco
import geopandas as gpd

from scripts._helpers import _simplify_polys, configure_logging, set_scenario_config

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()

EUROPE_COUNTRIES = [
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
    "UA",
    "MD",
]


def eez(eez_path, country_list=EUROPE_COUNTRIES):
    df = gpd.read_file(eez_path)
    iso3_list = cc.convert(country_list, src="ISO2", to="ISO3")  # noqa: F841
    pol_type = ["200NM", "Overlapping claim"]  # noqa: F841
    df = df.query("ISO_TER1 in @iso3_list and POL_TYPE in @pol_type").copy()
    df["name"] = cc.convert(df.ISO_TER1, src="ISO3", to="ISO2")
    s = df.set_index("name").geometry.map(
        lambda s: _simplify_polys(s, minarea=0.1, filterremote=False)
    )
    s = s.to_frame("geometry").set_crs(df.crs)
    s.index.name = "name"
    return s


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_offshore_shapes")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    offshore_shapes = eez(snakemake.input.eez, snakemake.params.countries)
    offshore_shapes.reset_index().to_file(snakemake.output.offshore_shapes)
