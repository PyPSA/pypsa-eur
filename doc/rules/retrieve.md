<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Retrieving Data {#data}

Not all data dependencies are shipped with the git repository, since git is not suited for handling large changing files.
Instead we use separate steps in the workflow (`rules` executed by `snakemake`) to download external data using the `retrieve_<dataset>` rules.

Data is generally retrieved in a version-controlled manner, enabling control over input data versions, reproducibility and consistency of modelling runs.
The rules download data into subfolders in the `data/` directory, following the structure
`data/{dataset}/{source}/{version}`, e.g. `data/jrc_idees/primary/March-2025-V1/`.
Which specific data version is retrieved can be controlled in the [data configuration](../configuration.md#data_cf).

Every dataset is listed with its owner, link and license in the
[data inventory](../data_sources.md#data-inventory). Its available versions and
sources (`archive`, `primary` or `build`) are registered in `data/versions.csv`.
Most `retrieve_<dataset>` rules simply download the registered URL for the
selected source. The rules below need more than that.

## Rules requiring credentials {#credentials}

Credentials are read from environment variables. You can put them into a `.env`
file in the repository root, which is loaded automatically and ignored by git.

| Rule | Requirement |
|------|-------------|
| [build_cutout][] | A [Copernicus Climate Data Store](https://cds.climate.copernicus.eu) account with the `cdsapi` key set up as described [on their website](https://cds.climate.copernicus.eu/how-to-api). Only needed when `data: cutout: source: build`. |
| `retrieve_electricity_demand_entsoe` | An [ENTSO-E Transparency Platform](https://transparency.entsoe.eu) API token in `ENTSOE_API_TOKEN`. Only needed when `data: entsoe_electricity_demand: source: build`. |
| `retrieve_corine` | A [Copernicus Land Monitoring Service](https://land.copernicus.eu/user/login) API key in `CORINE_API_TOKEN`. Only needed when `data: corine: source: primary`. |
| `retrieve_seawater_temperature` | [Copernicus Marine Service](https://marine.copernicus.eu/) credentials configured for the `copernicusmarine` package. Only needed when heat pumps use the `sea_water` heat source with a non-test cutout. |

## Rules with special handling {#special}

- `retrieve_cutout` downloads pre-built weather cutouts; see [cutouts](../configuration.md#atlite_cf).
- `retrieve_osm_data_raw` queries the [Overpass API](https://overpass-api.de) per country when `data: osm: source: build`. The endpoint, retries and user agent are set under [overpass_api](../configuration.md#overpass_api_cf).
- `retrieve_wdpa` and `retrieve_wdpa_marine` resolve the monthly changing download URL of the [World Database on Protected Areas](https://www.protectedplanet.net/). The data may not be redistributed, so the `archive` source points to a web archive copy.
- `retrieve_bidding_zones_entsoepy` and `retrieve_bidding_zones_electricitymaps` download bidding zone shapes via the [entsoe-py](https://github.com/EnergieID/entsoe-py) package and from [Electricity Maps](https://github.com/electricitymaps/electricitymaps-contrib). They are combined by [build_bidding_zones][].
- `retrieve_cost_data` downloads techno-economic assumptions from the [technology-data repository](https://github.com/pypsa/technology-data) as `data/costs/{source}/{version}/costs_{horizon}.csv`. The cost year can be fixed with `costs: year:` (see [costs](../configuration.md#costs_cf)).

## Electricity demand data {#demand}

Historical hourly electricity demand is retrieved from three sources and
combined by [build_electricity_demand][]:

| Rule | Source | Output |
|------|--------|--------|
| `retrieve_electricity_demand_opsd` | [OPSD platform](https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv), per country | `data/opsd_electricity_demand/{source}/{version}/electricity_demand_opsd_raw.csv` |
| `retrieve_electricity_demand_entsoe` | [ENTSO-E Transparency Platform](https://transparency.entsoe.eu), per country | `data/entsoe_electricity_demand/{source}/{version}/electricity_demand_entsoe_raw.csv` |
| `retrieve_electricity_demand_neso` | [NESO Data Portal](https://www.neso.energy/data-portal/historic-demand-data), United Kingdom | `data/neso_electricity_demand/{source}/{version}/electricity_demand_neso_raw.csv` |

The spatial distribution of demand within countries uses:

- `retrieve_electricity_demand_energy_atlas`: a 1 km by 1 km raster of estimated annual electricity demand from the [JRC Energy Atlas](https://energy-industry-geolab.jrc.ec.europa.eu/energy-atlas/).
- `retrieve_desnz_electricity_consumption`: subnational electricity consumption for Great Britain from the [Department for Energy Security and Net Zero](https://www.gov.uk/government/statistics/regional-and-local-authority-electricity-consumption-statistics).
- `retrieve_ons_lad`: shapefiles of local authorities in the United Kingdom from the [Office for National Statistics](https://geoportal.statistics.gov.uk/datasets/ons::local-authority-districts-may-2024-boundaries-uk-bsc-2/about).
