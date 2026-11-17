# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build the afforestation potential, sequestration rate and rotation age per region.

The NUTS2 sequestration rates and rotation ages are averaged per region,
weighted by the overlap area. The rates are means over the rotation, derived
from forest yield tables. Regions without NUTS2 data take the mean of their
country. The potential is the rate times the usable land of the region.

References
----------
- [Fernandes et al. (2026)](https://arxiv.org/abs/2603.25663), method and
  dataset.
- [Pilli et al. (2024)](https://publications.jrc.ec.europa.eu/repository/handle/JRC135639),
  forest growth data behind the sequestration rates.
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import area, configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

RENAME = {"CO2 seq rate tCO2/(ha y)": "rate", "rotation_age_years": "rotation_age"}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_afforestation_potentials")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    regions = gpd.read_file(snakemake.input.regions)
    land = pd.read_csv(snakemake.input.land, index_col=0)["area"]

    rename_countries = {"^EL": "GR", "^UK": "GB"}
    rates = pd.read_csv(snakemake.input.rates, index_col=0)[list(RENAME)]
    rates = rates.rename(columns=RENAME)
    rates.index = rates.index.to_series().replace(rename_countries, regex=True)

    nuts2 = gpd.read_file(snakemake.input.nuts2)
    nuts2["NUTS_ID"] = nuts2["NUTS_ID"].replace(rename_countries, regex=True)
    nuts2 = nuts2.set_index("NUTS_ID")[["geometry"]].join(rates, how="inner")

    overlay = gpd.overlay(regions[["name", "geometry"]], nuts2, keep_geom_type=True)
    weight = area(overlay)
    weighted = overlay[rates.columns].mul(weight, axis=0).groupby(overlay["name"]).sum()
    df = weighted.div(weight.groupby(overlay["name"]).sum(), axis=0)

    country_mean = rates.groupby(rates.index.str[:2]).mean()
    names = pd.Index(regions["name"])
    df = df.reindex(names).fillna(
        country_mean.reindex(names.str[:2]).set_axis(names),
    )

    # usable land in ha times rate in tCO2/ha/a
    df["potential"] = land * 100 * df["rate"]

    logger.info(f"Afforestation potential: {df['potential'].sum() / 1e6:.1f} MtCO2/a")
    df.to_csv(snakemake.output[0])
