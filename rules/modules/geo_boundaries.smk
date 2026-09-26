# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import country_converter as coco

GEO_BOUNDARIES_DIR = "resources/modules/geo_boundaries"


def geo_boundaries_scenario(cfg: dict) -> dict:
    """Translate the countries and module settings of one (scenario) config."""
    opts = cfg["modules"]["geo_boundaries"]
    nuts = {
        "source": "nuts",
        "subtype": "3",
        "resolution": opts["nuts_resolution"],
        "year": opts["nuts_year"],
    }
    adm1 = {"source": "geoboundaries", "subtype": "1", "release_type": "gbOpen"}
    countries = cfg["countries"]
    iso3 = coco.CountryConverter().pandas_convert(pd.Series(countries), to="ISO3")
    return {
        "countries": {
            a3: adm1 if a2 in opts["adm1_countries"] else nuts
            for a2, a3 in zip(countries, iso3)
        }
    }


if run["scenarios"]["enable"]:
    geo_boundaries_scenarios = {
        name: geo_boundaries_scenario(scenario_config(name)) for name in scenarios
    }
else:
    geo_boundaries_scenarios = {"default": geo_boundaries_scenario(config)}


module geo_boundaries:
    pathvars:
        logs="logs/modules/geo_boundaries",
        resources="data/modules/geo_boundaries",
        results=GEO_BOUNDARIES_DIR,
    snakefile:
        github(
            "modelblocks-org/module_geo_boundaries",
            path="workflow/Snakefile",
            tag=config["modules"]["geo_boundaries"]["version"],
        )
    config:
        {
            "crs": {"projected": "epsg:3035", "geographic": "epsg:4326"},
            "voronoi_eez": {"enabled": False},
            "scenarios": geo_boundaries_scenarios,
        }


use rule * from geo_boundaries exclude all as geo_boundaries_*


def geo_boundaries_shapes(w):
    """Module output for this run: one module scenario per PyPSA-Eur scenario."""
    return f"{GEO_BOUNDARIES_DIR}/{w.get('run', 'default')}/shapes.parquet"
