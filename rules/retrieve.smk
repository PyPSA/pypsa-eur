# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import os
import requests
from datetime import datetime
from dateutil.relativedelta import relativedelta
from shutil import move, unpack_archive, rmtree, copy2
from zipfile import ZipFile


# Configure the default storage provider for accessing remote files using http
# and the special storage plugin for accessing Zenodo files
storage:
    provider="http",
    keep_local=True,
    retries=3,


storage cached_http:
    provider="cached-http",


if (EUROSTAT_BALANCES_DATASET := dataset_version("eurostat_balances"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eurostat_balances:
        """Retrieves Eurostat energy balances by country and fuel."""
        input:
            tsv_gz=storage(EUROSTAT_BALANCES_DATASET["url"]),
        output:
            tsv_gz=f"{EUROSTAT_BALANCES_DATASET['folder']}/estat_nrg_bal_c.tsv.gz",
        run:
            copy2(input["tsv_gz"], output["tsv_gz"])


if (
    EUROSTAT_HOUSEHOLD_BALANCES_DATASET := dataset_version(
        "eurostat_household_balances"
    )
)["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eurostat_household_balances:
        """Retrieves Eurostat household final energy consumption balances."""
        input:
            csv=storage(EUROSTAT_HOUSEHOLD_BALANCES_DATASET["url"]),
        output:
            csv=f"{EUROSTAT_HOUSEHOLD_BALANCES_DATASET['folder']}/nrg_d_hhq.csv",
        run:
            copy2(input["csv"], output["csv"])


if (SWISS_ENERGY_BALANCES_DATASET := dataset_version("swiss_energy_balances"))[
    "source"
] in [
    "archive",
    "primary",
]:

    rule retrieve_swiss_energy_balances:
        """Retrieves Swiss energy balances from the Swiss Federal Office of Energy."""
        input:
            xlsx=storage(SWISS_ENERGY_BALANCES_DATASET["url"]),
        output:
            xlsx=f"{SWISS_ENERGY_BALANCES_DATASET['folder']}/12361-VWZ_Webtabellen_2024.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


if (NUTS3_POPULATION_DATASET := dataset_version("nuts3_population"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_nuts3_population:
        """Retrieves Eurostat population by NUTS3 region."""
        input:
            gz=storage(NUTS3_POPULATION_DATASET["url"]),
        output:
            gz=f"{NUTS3_POPULATION_DATASET['folder']}/nama_10r_3popgdp.tsv.gz",
        retries: 2
        run:
            copy2(input["gz"], output["gz"])


if (CORINE_DATASET := dataset_version("corine"))["source"] in ["archive"]:

    rule retrieve_corine:
        """Retrieves the CORINE land cover raster from the PyPSA data archive and unpacks it."""
        input:
            zip_file=storage(CORINE_DATASET["url"]),
        output:
            zip_file=f"{CORINE_DATASET['folder']}/corine.zip",
            tif_file=f"{CORINE_DATASET['folder']}/corine.tif",
        run:
            output_folder = Path(output["zip_file"]).parent
            unpack_archive(input["zip_file"], output_folder)
            copy2(input["zip_file"], output["zip_file"])
            copy2(
                f"{output_folder}/corine/g250_clc06_V18_5.tif", output["tif_file"]
            )

elif (CORINE_DATASET := dataset_version("corine"))["source"] in ["primary"]:

    rule retrieve_corine:
        """Downloads the CORINE land cover raster from the Copernicus Land Monitoring Service API."""
        output:
            zip=f"{CORINE_DATASET['folder']}/corine.zip",
            tif_file=f"{CORINE_DATASET['folder']}/corine.tif",
        log:
            logs("retrieve_corine_primary.log"),
        retries: 2
        resources:
            mem_mb=1000,
        params:
            apikey=os.environ.get("CORINE_API_TOKEN", ""),
        script:
            scripts("retrieve_corine_dataset_primary.py")


if (H2_SALT_CAVERNS_DATASET := dataset_version("h2_salt_caverns"))["source"] in [
    "archive"
]:

    rule retrieve_h2_salt_caverns:
        """Retrieves hydrogen salt cavern storage potentials in GWh per square kilometre."""
        input:
            geojson=storage(H2_SALT_CAVERNS_DATASET["url"]),
        output:
            geojson=f"{H2_SALT_CAVERNS_DATASET['folder']}/h2_salt_caverns_GWh_per_sqkm.geojson",
        retries: 2
        run:
            copy2(input["geojson"], output["geojson"])


if (GDP_PER_CAPITA_DATASET := dataset_version("gdp_per_capita"))["source"] in [
    "archive"
]:

    rule retrieve_gdp_per_capita:
        """Retrieves the gridded GDP per capita (PPP) dataset by Kummu et al."""
        input:
            gdp=storage(GDP_PER_CAPITA_DATASET["url"]),
        output:
            gdp=f"{GDP_PER_CAPITA_DATASET['folder']}/GDP_per_capita_PPP_1990_2015_v2.nc",
        retries: 2
        run:
            copy2(input["gdp"], output["gdp"])


if (POPULATION_COUNT_DATASET := dataset_version("population_count"))["source"] in [
    "archive",
    "primary",
]:

    rule retrieve_population_count:
        """Retrieves the WorldPop gridded population count raster."""
        input:
            tif=storage(POPULATION_COUNT_DATASET["url"]),
        output:
            tif=f"{POPULATION_COUNT_DATASET['folder']}/ppp_2019_1km_Aggregated.tif",
        retries: 2
        run:
            copy2(input["tif"], output["tif"])
            if POPULATION_COUNT_DATASET["source"] == "primary":
                import xarray as xr
                import rioxarray as rio

                file_path = output["tif"]
                ds = xr.open_dataarray(file_path)
                ds_reqd = ds.sel(x=slice(15.55, 40.41), y=slice(52.49, 41.72))
                ds_reqd.rio.to_raster(file_path)


if (GHG_EMISSIONS_DATASET := dataset_version("ghg_emissions"))["source"] in [
    "archive",
    "primary",
]:

    rule retrieve_ghg_emissions:
        """Retrieves national greenhouse gas emissions reported to the UNFCCC and the EEA."""
        input:
            ghg=storage(GHG_EMISSIONS_DATASET["url"]),
        output:
            csv=f"{GHG_EMISSIONS_DATASET['folder']}/UNFCCC_v23.csv",
            zip=(
                f"{GHG_EMISSIONS_DATASET['folder']}/UNFCCC_v23.csv.zip"
                if GHG_EMISSIONS_DATASET["source"] == "primary"
                else []
            ),
            directory=(
                directory(GHG_EMISSIONS_DATASET["folder"])
                if GHG_EMISSIONS_DATASET["source"] == "primary"
                else []
            ),
        retries: 2
        run:
            if GHG_EMISSIONS_DATASET["source"] == "primary":
                copy2(input["ghg"], output["zip"])
                unpack_archive(output["zip"], GHG_EMISSIONS_DATASET["folder"])
            else:
                copy2(input["ghg"], output["csv"])


if (GEBCO_DATASET := dataset_version("gebco"))["source"] in ["archive", "primary"]:

    rule retrieve_gebco:
        """Retrieves the GEBCO bathymetry grid used for offshore wind depth limits."""
        input:
            storage(GEBCO_DATASET["url"]),
        output:
            gebco=f"{GEBCO_DATASET['folder']}/GEBCO_2014_2D.nc",
            zip_file=(
                f"{GEBCO_DATASET['folder']}/GEBCO_2014.zip"
                if GEBCO_DATASET["source"] == "primary"
                else []
            ),
        run:
            if GEBCO_DATASET["source"] == "primary":
                import xarray as xr

                copy2(input[0], output["zip_file"])
                output_folder = Path(output["zip_file"]).parent
                unpack_archive(output["zip_file"], output_folder)
                # Limit extent to Europe to reduce file size
                ds = xr.open_dataset(output["gebco"])
                ds = ds.sel(lat=slice(32, 73), lon=slice(-21, 45))
                ds.to_netcdf(output["gebco"])
            else:
                copy2(input[0], output["gebco"])


if (ATTRIBUTED_PORTS_DATASET := dataset_version("attributed_ports"))["source"] in [
    "archive",
    "primary",
]:

    rule retrieve_attributed_ports:
        """Retrieves the World Bank dataset of international ports with attributes."""
        input:
            json=storage(ATTRIBUTED_PORTS_DATASET["url"]),
        output:
            json=f"{ATTRIBUTED_PORTS_DATASET['folder']}/attributed_ports.json",
        retries: 2
        run:
            copy2(input["json"], output["json"])


if (JRC_IDEES_DATASET := dataset_version("jrc_idees"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_jrc_idees:
        """Retrieves and unpacks the JRC-IDEES energy-economy-emissions dataset."""
        input:
            zip_file=storage(JRC_IDEES_DATASET["url"]),
        output:
            zip_file=f"{JRC_IDEES_DATASET['folder']}/jrc_idees.zip",
            directory=directory(JRC_IDEES_DATASET["folder"]),
        run:
            copy2(input["zip_file"], output["zip_file"])
            output_folder = Path(output["zip_file"]).parent
            unpack_archive(output["zip_file"], output_folder)


if (EU_NUTS2013_DATASET := dataset_version("eu_nuts2013"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eu_nuts_2013:
        """Retrieves and unpacks the Eurostat NUTS 2013 region shapes."""
        input:
            shapes=storage(EU_NUTS2013_DATASET["url"]),
        output:
            zip_file=f"{EU_NUTS2013_DATASET['folder']}/ref-nuts-2013-03m.geojson.zip",
            folder=directory(
                f"{EU_NUTS2013_DATASET['folder']}/ref-nuts-2013-03m.geojson"
            ),
            shapes_level_3=f"{EU_NUTS2013_DATASET['folder']}/ref-nuts-2013-03m.geojson/NUTS_RG_03M_2013_4326_LEVL_3.geojson",
            shapes_level_2=f"{EU_NUTS2013_DATASET['folder']}/ref-nuts-2013-03m.geojson/NUTS_RG_03M_2013_4326_LEVL_2.geojson",
        run:
            copy2(input["shapes"], output["zip_file"])
            unpack_archive(output["zip_file"], Path(output.shapes_level_3).parent)


if (EU_NUTS2021_DATASET := dataset_version("eu_nuts2021"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eu_nuts_2021:
        """Retrieves and unpacks the Eurostat NUTS 2021 region shapes."""
        input:
            shapes=storage(EU_NUTS2021_DATASET["url"]),
        output:
            zip_file=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson.zip",
            folder=directory(
                f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson"
            ),
            shapes_level_3=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_3.geojson",
            shapes_level_2=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_2.geojson",
            shapes_level_1=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_1.geojson",
            shapes_level_0=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_0.geojson",
        run:
            copy2(input["shapes"], output["zip_file"])
            unpack_archive(output["zip_file"], Path(output.shapes_level_3).parent)


if (
    BIDDING_ZONES_ELECTRICITYMAPS_DATASET := dataset_version(
        "bidding_zones_electricitymaps"
    )
)["source"] in ["primary", "archive"]:

    rule retrieve_bidding_zones_electricitymaps:
        """Retrieves bidding zone shapes from Electricity Maps."""
        input:
            geojson=storage(BIDDING_ZONES_ELECTRICITYMAPS_DATASET["url"]),
        output:
            geojson=f"{BIDDING_ZONES_ELECTRICITYMAPS_DATASET['folder']}/bidding_zones_electricitymaps.geojson",
        log:
            "logs/retrieve_bidding_zones_electricitymaps.log",
        retries: 2
        resources:
            mem_mb=1000,
        run:
            copy2(input["geojson"], output["geojson"])


if (BIDDING_ZONES_ENTSOEPY_DATASET := dataset_version("bidding_zones_entsoepy"))[
    "source"
] in ["primary", "archive"]:

    rule retrieve_bidding_zones_entsoepy:
        """Downloads bidding zone shapes for all entsoe-py areas and merges them into one file."""
        output:
            geojson=f"{BIDDING_ZONES_ENTSOEPY_DATASET['folder']}/bidding_zones_entsoepy.geojson",
        log:
            "logs/retrieve_bidding_zones_entsoepy.log",
        retries: 2
        resources:
            mem_mb=1000,
        run:
            import entsoe
            import geopandas as gpd
            from urllib.error import HTTPError, URLError

            logger.info("Downloading entsoe-py zones...")
            gdfs: list[gpd.GeoDataFrame] = []
            url = f"{BIDDING_ZONES_ENTSOEPY_DATASET['url']}"
            for area in entsoe.Area:
                name = area.name
                try:
                    file_url = f"{url}/{name}.geojson"
                    gdfs.append(gpd.read_file(file_url))
                except HTTPError as e:
                    logger.debug(f"Area file not available for {name}: {e}")
                    continue
                except (URLError, TimeoutError) as e:
                    raise Exception(f"Network error retrieving {name}: {e}")
            shapes = pd.concat(gdfs, ignore_index=True)  # type: ignore
            logger.info("Downloading entsoe-py zones... Done")
            shapes.to_file(output.geojson)


if (CUTOUT_DATASET := dataset_version("cutout"))["source"] in [
    "archive",
]:

    rule retrieve_cutout:
        """Retrieves pre-built atlite weather cutouts from the PyPSA data archive."""
        input:
            storage(CUTOUT_DATASET["url"] + "/{cutout}.nc"),
        output:
            CUTOUT_DATASET["folder"] + "/{cutout}.nc",
        log:
            "logs/retrieve_cutout/{cutout}.log",
        retries: 2
        resources:
            mem_mb=5000,
        run:
            copy2(input[0], output[0])


if (COUNTRY_RUNOFF_DATASET := dataset_version("country_runoff"))["source"] in [
    "archive"
]:

    rule retrieve_country_runoff:
        """Retrieves country-level daily runoff sums derived from ERA5."""
        input:
            storage(COUNTRY_RUNOFF_DATASET["url"]),
        output:
            era5_runoff=f"{COUNTRY_RUNOFF_DATASET['folder']}/era5-runoff-per-country.csv",
        run:
            copy2(input[0], output[0])


if (COUNTRY_HDD_DATASET := dataset_version("country_hdd"))["source"] in ["archive"]:

    rule retrieve_country_hdd:
        """Retrieves country-level heating degree days derived from ERA5."""
        input:
            storage(COUNTRY_HDD_DATASET["url"]),
        output:
            era5_runoff=f"{COUNTRY_HDD_DATASET['folder']}/era5-HDD-per-country.csv",
        run:
            copy2(input[0], output[0])


if (COSTS_DATASET := dataset_version("costs"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_cost_data:
        """Retrieves technology cost assumptions for a planning horizon from technology-data."""
        input:
            costs=storage(COSTS_DATASET["url"] + "/costs_{horizon}.csv"),
        output:
            costs=COSTS_DATASET["folder"] + "/costs_{horizon}.csv",
        run:
            copy2(input["costs"], output["costs"])


if (POWERPLANTS_DATASET := dataset_version("powerplants"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_powerplants:
        """Retrieves the powerplantmatching dataset of European power plants."""
        input:
            powerplants=storage(POWERPLANTS_DATASET["url"]),
        output:
            powerplants=f"{POWERPLANTS_DATASET['folder']}/powerplants.csv",
        run:
            copy2(input["powerplants"], output["powerplants"])


if (SCIGRID_GAS_DATASET := dataset_version("scigrid_gas"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_gas_infrastructure_data:
        """Retrieves and unpacks the SciGRID_gas IGGIELGN gas transmission network."""
        input:
            zip_file=storage(SCIGRID_GAS_DATASET["url"]),
        output:
            zip_file=f"{SCIGRID_GAS_DATASET['folder']}/IGGIELGN.zip",
            entry=f"{SCIGRID_GAS_DATASET['folder']}/data/IGGIELGN_BorderPoints.geojson",
            storage=f"{SCIGRID_GAS_DATASET['folder']}/data/IGGIELGN_Storages.geojson",
            gas_network=f"{SCIGRID_GAS_DATASET['folder']}/data/IGGIELGN_PipeSegments.geojson",
        run:
            copy2(input["zip_file"], output["zip_file"])
            output_folder = Path(output["zip_file"]).parent
            unpack_archive(output["zip_file"], output_folder)


if (OPSD_DEMAND_DATA := dataset_version("opsd_electricity_demand"))["source"] in [
    "build"
]:

    rule retrieve_electricity_demand_opsd:
        """Builds the OPSD electricity demand time series from the OPSD data platform."""
        output:
            csv=f"{OPSD_DEMAND_DATA['folder']}/electricity_demand_opsd_raw.csv",
        log:
            "logs/retrieve_electricity_demand_opsd.log",
        retries: 2
        resources:
            mem_mb=5000,
        params:
            versions=["2019-06-05", "2020-10-06"],
        script:
            scripts("retrieve_electricity_demand_opsd.py")


if (OPSD_DEMAND_DATA := dataset_version("opsd_electricity_demand"))["source"] in [
    "archive"
]:

    rule retrieve_electricity_demand_opsd:
        """Retrieves the OPSD electricity demand time series from the PyPSA data archive."""
        input:
            csv=storage(OPSD_DEMAND_DATA["url"]),
        output:
            csv=f"{OPSD_DEMAND_DATA['folder']}/electricity_demand_opsd_raw.csv",
        retries: 2
        run:
            copy2(input["csv"], output["csv"])


if (ENTSOE_DEMAND_DATA := dataset_version("entsoe_electricity_demand"))["source"] in [
    "build"
]:

    ENTSOE_COUNTRIES = [
        "AL",
        "AT",
        "BE",
        "BA",
        "BG",
        "CH",
        "CY",
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
        "MD",
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
        "UA",
        "XK",
    ]

    rule retrieve_electricity_demand_entsoe_country:
        """Downloads electricity demand time series for one country from the ENTSO-E Transparency Platform."""
        output:
            csv=f"{ENTSOE_DEMAND_DATA['folder']}"
            + "/electricity_demand_entsoe_raw_{country}.csv",
        log:
            "logs/retrieve_electricity_demand_entsoe_{country}.log",
        retries: 2
        resources:
            mem_mb=2000,
        params:
            entsoe_token=os.environ.get("ENTSOE_API_TOKEN", ""),
        script:
            scripts("retrieve_electricity_demand_entsoe.py")

    rule retrieve_electricity_demand_entsoe:
        """Merges per-country ENTSO-E electricity demand time series into one file."""
        input:
            csvs=expand(
                f"{ENTSOE_DEMAND_DATA['folder']}"
                + "/electricity_demand_entsoe_raw_{country}.csv",
                country=ENTSOE_COUNTRIES,
            ),
        output:
            csv=f"{ENTSOE_DEMAND_DATA['folder']}/electricity_demand_entsoe_raw.csv",
        run:
            import pandas as pd

            loads = [pd.read_csv(csv, index_col=0) for csv in input.csvs]
            df = pd.concat(loads, axis=1, join="outer").sort_index()
            df.to_csv(output.csv)


if (ENTSOE_DEMAND_DATA := dataset_version("entsoe_electricity_demand"))["source"] in [
    "archive"
]:

    rule retrieve_electricity_demand_entsoe:
        """Retrieves the ENTSO-E electricity demand time series from the PyPSA data archive."""
        input:
            csv=storage(ENTSOE_DEMAND_DATA["url"]),
        output:
            csv=f"{ENTSOE_DEMAND_DATA['folder']}/electricity_demand_entsoe_raw.csv",
        retries: 2
        run:
            copy2(input["csv"], output["csv"])


if (NESO_DEMAND_DATA := dataset_version("neso_electricity_demand"))["source"] in [
    "build"
]:

    rule retrieve_electricity_demand_neso:
        """Downloads Great Britain electricity demand time series from the NESO data portal."""
        output:
            csv=f"{NESO_DEMAND_DATA['folder']}/electricity_demand_neso_raw.csv",
        log:
            "logs/retrieve_electricity_demand_neso.log",
        retries: 2
        resources:
            mem_mb=5000,
        script:
            scripts("retrieve_electricity_demand_neso.py")


if (NESO_DEMAND_DATA := dataset_version("neso_electricity_demand"))["source"] in [
    "archive"
]:

    rule retrieve_electricity_demand_neso:
        """Retrieves the NESO electricity demand time series from the PyPSA data archive."""
        input:
            csv=storage(NESO_DEMAND_DATA["url"]),
        output:
            csv=f"{NESO_DEMAND_DATA['folder']}/electricity_demand_neso_raw.csv",
        retries: 2
        run:
            copy2(input["csv"], output["csv"])


if (
    SYNTHETIC_ELECTRICITY_DEMAND_DATASET := dataset_version(
        "synthetic_electricity_demand"
    )
)["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_synthetic_electricity_demand:
        """Retrieves synthetic hourly electricity demand time series by country."""
        input:
            csv=storage(SYNTHETIC_ELECTRICITY_DEMAND_DATASET["url"]),
        output:
            csv=f"{SYNTHETIC_ELECTRICITY_DEMAND_DATASET['folder']}/load_synthetic_raw.csv",
        retries: 2
        run:
            copy2(input["csv"], output["csv"])


if (ENERGY_ATLAS_DATASET := dataset_version("jrc_energy_atlas"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_electricity_demand_energy_atlas:
        """Downloads the JRC Energy Atlas raster of annual electricity demand."""
        output:
            tif=f"{ENERGY_ATLAS_DATASET['folder']}/electricity_tot_demand_2019.tif",
        run:
            import requests

            url = ENERGY_ATLAS_DATASET["url"]
            response = requests.get(url)
            response.raise_for_status()
            with open(output["tif"], "wb") as f:
                f.write(response.content)


if (
    DESNZ_ELECTRICITY_CONSUMPTION_DATASET := dataset_version(
        "desnz_electricity_consumption"
    )
)["source"] in ["primary", "archive"]:

    rule retrieve_desnz_electricity_consumption:
        """Downloads DESNZ subnational electricity consumption statistics for the UK."""
        output:
            xlsx=f"{DESNZ_ELECTRICITY_CONSUMPTION_DATASET['folder']}/Subnational_electricity_consumption_statistics_2005-2024.xlsx",
        run:
            import requests

            url = DESNZ_ELECTRICITY_CONSUMPTION_DATASET["url"]
            response = requests.get(url)
            response.raise_for_status()
            with open(output["xlsx"], "wb") as f:
                f.write(response.content)


if (ONS_LAD_DATASET := dataset_version("ons_lad"))["source"] in ["archive"]:

    rule retrieve_ons_lad:
        """Retrieves UK Local Authority District boundaries from the PyPSA data archive."""
        input:
            geojson=storage(ONS_LAD_DATASET["url"]),
        output:
            geojson=f"{ONS_LAD_DATASET['folder']}/Local_Authority_Districts_May_2024_Boundaries__UK_BSC.geojson",
        run:
            copy2(input["geojson"], output["geojson"])

elif ONS_LAD_DATASET["source"] in ["primary"]:

    rule retrieve_ons_lad:
        """Downloads UK Local Authority District boundaries from the ONS ArcGIS service."""
        output:
            geojson=f"{ONS_LAD_DATASET['folder']}/Local_Authority_Districts_May_2024_Boundaries__UK_BSC.geojson",
        run:
            import requests

            url = ONS_LAD_DATASET["url"]
            params = {
                "outFields": "*",
                "where": "1=1",
                "f": "geojson",
            }
            response = requests.get(url, params=params)
            with open(output["geojson"], "wb") as f:
                f.write(response.content)


if (SHIP_RASTER_DATASET := dataset_version("ship_raster"))["source"] in [
    "archive",
    "primary",
]:

    rule retrieve_ship_raster:
        """Retrieves the World Bank global shipping traffic density raster."""
        input:
            zip_file=storage(SHIP_RASTER_DATASET["url"]),
        output:
            zip_file=f"{SHIP_RASTER_DATASET['folder']}/shipdensity_global.zip",
        log:
            "logs/retrieve_ship_raster.log",
        retries: 2
        resources:
            mem_mb=5000,
        run:
            copy2(input["zip_file"], output["zip_file"])


if (ENSPRESO_BIOMASS_DATASET := dataset_version("enspreso_biomass"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_enspreso_biomass:
        """Retrieves JRC ENSPRESO biomass potentials."""
        input:
            xlsx=storage(ENSPRESO_BIOMASS_DATASET["url"]),
        output:
            xlsx=f"{ENSPRESO_BIOMASS_DATASET['folder']}/ENSPRESO_BIOMASS.xlsx",
        retries: 1
        run:
            copy2(input["xlsx"], output["xlsx"])


if (TABULA_CALCULATOR := dataset_version("tabula_calculator"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_tabula_calculator:
        """Retrieves the TABULA building typology calculator workbook."""
        input:
            xlsx=storage(TABULA_CALCULATOR["url"]),
        output:
            xlsx=f"{TABULA_CALCULATOR['folder']}/tabula-calculator.xlsx",
        retries: 2
        run:
            copy2(input["xlsx"], output["xlsx"])


if (HOTMAPS_INDUSTRIAL_SITES := dataset_version("hotmaps_industrial_sites"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_hotmaps_industrial_sites:
        """Retrieves the Hotmaps database of industrial sites."""
        input:
            csv=storage(HOTMAPS_INDUSTRIAL_SITES["url"]),
        output:
            csv=f"{HOTMAPS_INDUSTRIAL_SITES['folder']}/Industrial_Database.csv",
        retries: 1
        run:
            copy2(input["csv"], output["csv"])


if (NITROGEN_STATISTICS_DATASET := dataset_version("nitrogen_statistics"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_nitrogen_statistics:
        """Retrieves USGS nitrogen supply and demand statistics."""
        input:
            xlsx=storage(NITROGEN_STATISTICS_DATASET["url"]),
        output:
            xlsx=f"{NITROGEN_STATISTICS_DATASET['folder']}/nitro-ert.xlsx",
        retries: 1
        run:
            copy2(input["xlsx"], output["xlsx"])


if (COPERNICUS_LAND_COVER_DATASET := dataset_version("copernicus_land_cover"))[
    "source"
] in ["primary", "archive"]:

    # Downloading Copernicus Global Land Cover for land cover and land use:
    # Website: https://land.copernicus.eu/global/products/lc
    rule download_copernicus_land_cover:
        """Retrieves the Copernicus Global Land Cover raster."""
        input:
            tif=storage(COPERNICUS_LAND_COVER_DATASET["url"]),
        output:
            tif=f"{COPERNICUS_LAND_COVER_DATASET['folder']}/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        run:
            copy2(input["tif"], output["tif"])


if (LUISA_LAND_COVER_DATASET := dataset_version("luisa_land_cover"))["source"] in [
    "primary",
    "archive",
]:

    # Downloading LUISA Base Map for land cover and land use:
    # Website: https://ec.europa.eu/jrc/en/luisa
    rule retrieve_luisa_land_cover:
        """Retrieves the JRC LUISA base map land cover raster."""
        input:
            tif=storage(LUISA_LAND_COVER_DATASET["url"]),
        output:
            tif=f"{LUISA_LAND_COVER_DATASET['folder']}/LUISA_basemap_020321_50m.tif",
        run:
            copy2(input["tif"], output["tif"])


if (EEZ_DATASET := dataset_version("eez"))["source"] in ["primary"]:

    rule retrieve_eez:
        """Downloads and unpacks the Marine Regions World EEZ shapes via the registration form."""
        output:
            zip_file=f"{EEZ_DATASET['folder']}/World_EEZ_{EEZ_DATASET['version']}_LR.zip",
            gpkg=f"{EEZ_DATASET['folder']}/World_EEZ_{EEZ_DATASET['version']}_LR/eez_{EEZ_DATASET['version'].split('_')[0]}_lowres.gpkg",
        run:
            from uuid import uuid4

            name = str(uuid4())[:8]
            org = str(uuid4())[:8]
            response = requests.post(
                f"{EEZ_DATASET['url']}",
                params={"name": f"World_EEZ_{EEZ_DATASET['version']}_LR.zip"},
                data={
                    "name": name,
                    "organisation": org,
                    "email": f"{name}@{org}.org",
                    "country": "Germany",
                    "user_category": "academia",
                    "purpose_category": "Research",
                    "agree": "1",
                },
            )
            with open(output["zip_file"], "wb") as f:
                f.write(response.content)
            output_folder = Path(output["zip_file"]).parent
            unpack_archive(output["zip_file"], output_folder)

elif (EEZ_DATASET := dataset_version("eez"))["source"] in ["archive"]:

    rule retrieve_eez:
        """Retrieves and unpacks the Marine Regions World EEZ shapes from the PyPSA data archive."""
        input:
            zip_file=storage(
                EEZ_DATASET["url"],
            ),
        output:
            zip_file=f"{EEZ_DATASET['folder']}/World_EEZ_{EEZ_DATASET['version']}_LR.zip",
            gpkg=f"{EEZ_DATASET['folder']}/World_EEZ_{EEZ_DATASET['version']}_LR/eez_{EEZ_DATASET['version'].split('_')[0]}_lowres.gpkg",
        run:
            output_folder = Path(output["zip_file"]).parent
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)


if (WB_URB_POP_DATASET := dataset_version("worldbank_urban_population"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_worldbank_urban_population:
        """Retrieves and unpacks the World Bank urban population share by country."""
        input:
            zip=storage(WB_URB_POP_DATASET["url"]),
        output:
            zip=f"{WB_URB_POP_DATASET['folder']}/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.zip",
            csv=f"{WB_URB_POP_DATASET['folder']}/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.csv",
        run:
            copy2(input["zip"], output["zip"])
            unpack_archive(output["zip"], WB_URB_POP_DATASET["folder"])
            # Filename contains some added numbers when downloaded,
            # remove them to have a consistent filename across versions
            target_filename = Path(output["csv"])
            origin_filename = next(
                Path(WB_URB_POP_DATASET["folder"]).rglob(
                    target_filename.stem + "*" + target_filename.suffix
                )
            )
            origin_filename.rename(output.csv)


if (CO2STOP_DATASET := dataset_version("co2stop"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_co2stop:
        """Retrieves and unpacks the JRC CO2Stop CO2 storage potentials."""
        input:
            zip_file=storage(CO2STOP_DATASET["url"]),
        output:
            zip_file=f"{CO2STOP_DATASET['folder']}/co2jrc_openformats.zip",
            storage_table=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Storage_Units.csv",
            storage_map=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_Polygons Data/StorageUnits_March13.kml",
            traps_table1=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps.csv",
            traps_table2=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps_Temp.csv",
            traps_table3=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps1.csv",
            traps_map=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_Polygons Data/DaughterUnits_March13.kml",
        run:
            output_folder = Path(output["zip_file"]).parent
            output_folder.mkdir(parents=True, exist_ok=True)
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)


if (GEM_EUROPE_GAS_TRACKER_DATASET := dataset_version("gem_europe_gas_tracker"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_gem_europe_gas_tracker:
        """Retrieves the Global Energy Monitor Europe Gas Tracker."""
        input:
            xlsx=storage(GEM_EUROPE_GAS_TRACKER_DATASET["url"]),
        output:
            xlsx="data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


if (GEM_GSPT_DATASET := dataset_version("gem_gspt"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_gem_steel_plant_tracker:
        """Retrieves the Global Energy Monitor Global Steel Plant Tracker."""
        input:
            xlsx=storage(GEM_GSPT_DATASET["url"]),
        output:
            xlsx=f"{GEM_GSPT_DATASET['folder']}/Global-Steel-Plant-Tracker.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


if (GEM_GCCT_DATASET := dataset_version("gem_gcct"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_gem_cement_concrete_tracker:
        """Retrieves the Global Energy Monitor Global Cement and Concrete Tracker."""
        input:
            xlsx=storage(GEM_GCCT_DATASET["url"]),
        output:
            xlsx=f"{GEM_GCCT_DATASET['folder']}/Global-Cement-and-Concrete-Tracker.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


if (BFS_ROAD_VEHICLE_STOCK_DATASET := dataset_version("bfs_road_vehicle_stock"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_bfs_road_vehicle_stock:
        """Retrieves the Swiss road vehicle stock from the Swiss Federal Statistical Office."""
        input:
            csv=storage(BFS_ROAD_VEHICLE_STOCK_DATASET["url"]),
        output:
            csv=f"{BFS_ROAD_VEHICLE_STOCK_DATASET['folder']}/vehicle_stock.csv",
        run:
            copy2(input["csv"], output["csv"])


if (BFS_GDP_AND_POPULATION_DATASET := dataset_version("bfs_gdp_and_population"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_bfs_gdp_and_population:
        """Retrieves Swiss GDP and population data from the Swiss Federal Statistical Office."""
        input:
            xlsx=storage(BFS_GDP_AND_POPULATION_DATASET["url"]),
        output:
            xlsx=f"{BFS_GDP_AND_POPULATION_DATASET['folder']}/gdp_and_population.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


def get_wdpa_url(DATASET) -> str:
    """
    Find the right URL for the WDPA / WDPA marine dataset based on the source type.
    """
    if DATASET["source"] == "archive":
        return DATASET["url"]
    elif DATASET["source"] == "primary":
        # Some logic to find the correct file URL from the WDPA website (primary source)
        # Sometimes files are released delayed or ahead of schedule, check which file is currently available
        def check_file_exists(url):
            response = requests.head(url)
            return response.status_code == 200

        # Basic pattern where WDPA files can be found
        url_pattern = DATASET["url"]

        # 3-letter month + 4 digit year for current/previous/next/pprevious/nnext months to test
        # order reflects priority of testing
        months = [
            datetime.now(),  # current
            (datetime.now() + relativedelta(months=-1)),  # previous month
            (datetime.now() + relativedelta(months=+1)),  # next month
            (datetime.now() + relativedelta(months=-2)),  # two months ago
            (datetime.now() + relativedelta(months=+2)),  # two months ahead
        ]
        months = [m.strftime("%b%Y") for m in months]

        # Test prioritised: current month -> previous -> next
        for bYYYY in months:
            url = url_pattern.format(bYYYY=bYYYY)
            if check_file_exists(url):
                return url

        raise ValueError(
            f"No {DATASET.dataset} files found at {url_pattern} for bY={months}."
        )


if (WDPA_DATASET := dataset_version("wdpa"))["source"] in [
    "primary",
    "archive",
]:

    # Downloading protected area database from WDPA
    # extract the main zip and then merge the contained 3 zipped shapefiles
    # Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
    rule retrieve_wdpa:
        """Retrieves the WDPA protected areas and merges the zipped shapefiles into one geopackage."""
        input:
            zip_file=storage(get_wdpa_url(WDPA_DATASET)),
        output:
            zip_file=f"{WDPA_DATASET['folder']}/WDPA_shp.zip",
            gpkg=f"{WDPA_DATASET['folder']}/WDPA.gpkg",
        retries: 2
        run:
            output_folder = Path(output["zip_file"]).parent
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)
            # Extract {bYYYY} from the input file / URL
            bYYYY = re.search(
                r"WDPA_(\w{3}\d{4})_Public_shp.zip",
                input["zip_file"],
            ).group(1)
            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = (
                    f"/vsizip/{output_folder}/WDPA_{bYYYY}_Public_shp_{i}.zip"
                )
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")


if (WDPA_MARINE_DATASET := dataset_version("wdpa_marine"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_wdpa_marine:
        # Downloading Marine protected area database from WDPA
        # extract the main zip and then merge the contained 3 zipped shapefiles
        # Website: https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
        """Retrieves the WDPA marine protected areas and merges the zipped shapefiles into one geopackage."""
        input:
            zip_file=storage(get_wdpa_url(WDPA_MARINE_DATASET)),
        output:
            zip_file=f"{WDPA_MARINE_DATASET['folder']}/WDPA_WDOECM_marine.zip",
            gpkg=f"{WDPA_MARINE_DATASET['folder']}/WDPA_WDOECM_marine.gpkg",
        retries: 2
        run:
            output_folder = Path(output["zip_file"]).parent
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)
            # Extract {bYYYY} from the input file / URL
            bYYYY = re.search(
                r"WDPA_WDOECM_(\w{3}\d{4})_Public_marine_shp.zip",
                input["zip_file"],
            ).group(1)
            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = f"/vsizip/{output_folder}/WDPA_WDOECM_{bYYYY}_Public_marine_shp_{i}.zip"
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")


if (INSTRAT_CO2_PRICES_DATASET := dataset_version("instrat_co2_prices"))["source"] in [
    "primary",
]:

    rule retrieve_co2_prices:
        """Downloads EU ETS CO2 allowance price time series from the Instrat energy API."""
        output:
            csv=f"{INSTRAT_CO2_PRICES_DATASET['folder']}/prices_eu_ets_all.csv",
        log:
            "logs/retrieve_co2_prices.log",
        retries: 2
        resources:
            mem_mb=5000,
        run:
            from io import StringIO
            import pandas as pd

            url = "https://energy-api.instrat.pl/api/prices/co2?all=1"
            headers = {
                "User-Agent": "Mozilla/5.0",
                "Accept": "application/json",
                "Referer": "https://energy.instrat.pl/",
            }
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            df = pd.read_json(StringIO(r.text))
            df.to_csv(output["csv"], index=False)


if (
    WORLD_BANK_COMMODITY_PRICES_DATASET := dataset_version("worldbank_commodity_prices")
)["source"] in ["primary", "archive"]:

    rule retrieve_worldbank_commodity_prices:
        """Retrieves the World Bank monthly commodity price time series."""
        input:
            xlsx=storage(WORLD_BANK_COMMODITY_PRICES_DATASET["url"]),
        output:
            xlsx=f"{WORLD_BANK_COMMODITY_PRICES_DATASET['folder']}/CMO-Historical-Data-Monthly.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


if (TYNDP_DATASET := dataset_version("tyndp"))["source"] in ["primary", "archive"]:

    rule retrieve_tyndp:
        """Retrieves and unpacks the ENTSO-E TYNDP reference grid and node lists."""
        input:
            line_data=storage(TYNDP_DATASET["url"] + "/Line-data.zip"),
            nodes=storage(TYNDP_DATASET["url"] + "/Nodes.zip"),
        output:
            line_data_zip=f"{TYNDP_DATASET['folder']}/Line-data.zip",
            nodes_zip=f"{TYNDP_DATASET['folder']}/Nodes.zip",
            reference_grid=f"{TYNDP_DATASET['folder']}/Line data/ReferenceGrid_Electricity.xlsx",
            nodes=f"{TYNDP_DATASET['folder']}/Nodes/LIST OF NODES.xlsx",
        log:
            "logs/retrieve_tyndp.log",
        run:
            for key in input.keys():
                # Keep zip file
                copy2(input[key], output[f"{key}_zip"])
                # unzip
                output_folder = Path(output[f"{key}_zip"]).parent
                unpack_archive(output[f"{key}_zip"], output_folder)
                # Remove __MACOSX directory if it exists
                macosx_dir = output_folder / "__MACOSX"
                rmtree(macosx_dir, ignore_errors=True)


def get_osm_archive_files(version):
    return [
        "buses.csv",
        "converters.csv",
        "lines.csv",
        "links.csv",
        "transformers.csv",
        # Newer versions include the additional map.html file for visualisation
        *(["map.html"] if float(version) >= 0.6 else []),
    ]


def input_base_network_incumbent(w):
    version = config_provider("osm_network_release", "compare_to", "version")(w)
    source = config_provider("osm_network_release", "compare_to", "source")(w)
    osm_dataset = dataset_version("osm", version=version, source=source)
    osm_path = osm_dataset["folder"]
    components = {"buses", "lines", "links", "converters", "transformers"}
    inputs = {c: f"{osm_path}/{c}.csv" for c in components}
    return inputs


if OSM_DATASET["source"] in ["archive"]:
    OSM_ARCHIVE_FILES = get_osm_archive_files(OSM_DATASET["version"])

    rule retrieve_osm_archive:
        """Retrieves the prebuilt OSM transmission grid from the PyPSA data archive."""
        input:
            **{
                file: storage(f"{OSM_DATASET['url']}/{file}")
                for file in OSM_ARCHIVE_FILES
            },
        output:
            **{file: f"{OSM_DATASET['folder']}/{file}" for file in OSM_ARCHIVE_FILES},
        log:
            "logs/retrieve_osm_archive.log",
        threads: 1
        resources:
            mem_mb=500,
        run:
            for key in input.keys():
                copy2(input[key], output[key])


# Only create incumbent rule if it points to a different folder
OSM_DATASET_INCUMBENT = dataset_version(
    "osm",
    version=config.get("osm_network_release", {})
    .get("compare_to", {})
    .get("version", "latest"),
    source=config.get("osm_network_release", {})
    .get("compare_to", {})
    .get("source", "archive"),
)

if OSM_DATASET_INCUMBENT["source"] in ["archive"] and OSM_DATASET_INCUMBENT[
    "folder"
] != OSM_DATASET.get("folder"):

    OSM_ARCHIVE_FILES_INCUMBENT = get_osm_archive_files(
        OSM_DATASET_INCUMBENT["version"]
    )

    rule retrieve_osm_archive_incumbent:
        """Retrieves a second OSM transmission grid release used for comparison with the current one."""
        input:
            **{
                file: storage(f"{OSM_DATASET_INCUMBENT['url']}/{file}")
                for file in OSM_ARCHIVE_FILES_INCUMBENT
            },
        output:
            **{
                file: f"{OSM_DATASET_INCUMBENT['folder']}/{file}"
                for file in OSM_ARCHIVE_FILES_INCUMBENT
            },
        log:
            "logs/retrieve_osm_archive_incumbent.log",
        threads: 1
        resources:
            mem_mb=500,
        run:
            for key in input.keys():
                copy2(input[key], output[key])


if OSM_DATASET["source"] == "build":
    OSM_RAW_JSON = [
        "cables_way.json",
        "lines_way.json",
        "routes_relation.json",
        "substations_way.json",
        "substations_relation.json",
    ]

    rule retrieve_osm_data_raw:
        """Downloads raw OSM power grid elements for one country via the Overpass API."""
        output:
            **{
                file.replace(
                    ".json", ""
                ): f"{OSM_DATASET['folder']}/{{country}}/{file}"
                for file in OSM_RAW_JSON
            },
        log:
            "logs/retrieve_osm_data_{country}.log",
        threads: 1
        params:
            overpass_api=config_provider("overpass_api"),
        script:
            scripts("retrieve_osm_data.py")

    rule retrieve_osm_data_raw_all:
        """Collects the raw OSM power grid data for all configured countries."""
        input:
            expand(
                f"{OSM_DATASET['folder']}/{{country}}/{{file}}",
                country=config_provider("countries"),
                file=OSM_RAW_JSON,
            ),


if (NATURA_DATASET := dataset_version("natura"))["source"] in ["archive"]:

    rule retrieve_natura:
        """Retrieves the prebuilt Natura 2000 raster."""
        input:
            storage(NATURA_DATASET["url"]),
        output:
            f"{NATURA_DATASET['folder']}/natura.tiff",
        log:
            "logs/retrieve_natura.log",
        run:
            copy2(input[0], output[0])

elif NATURA_DATASET["source"] == "build":

    rule build_natura_raster:
        """Downloads the Natura 2000 shapes from the EEA and rasterises them onto the cutout grid."""
        input:
            online=storage(NATURA_DATASET["url"]),
            cutout=lambda w: input_cutout(w),
        output:
            zip=f"{NATURA_DATASET['folder']}/raw/natura.zip",
            raw=directory(f"{NATURA_DATASET['folder']}/raw"),
            raster=f"{NATURA_DATASET['folder']}/natura.tiff",
        log:
            "logs/build_natura.log",
        resources:
            mem_mb=5000,
        script:
            scripts("build_natura.py")


if (OSM_BOUNDARIES_DATASET := dataset_version("osm_boundaries"))["source"] in [
    "primary"
]:

    rule retrieve_osm_boundaries:
        """Downloads OSM administrative boundaries for one country via the Overpass API."""
        output:
            json=f"{OSM_BOUNDARIES_DATASET['folder']}/{country}_adm1.json",
        log:
            "logs/retrieve_osm_boundaries_{country}_adm1.log",
        threads: 1
        script:
            scripts("retrieve_osm_boundaries.py")

elif (OSM_BOUNDARIES_DATASET := dataset_version("osm_boundaries"))["source"] in [
    "archive"
]:

    rule retrieve_osm_boundaries:
        """Retrieves and unpacks OSM administrative boundaries from the PyPSA data archive."""
        input:
            storage(
                f"{OSM_BOUNDARIES_DATASET['url']}",
            ),
        output:
            json1=f"{OSM_BOUNDARIES_DATASET['folder']}/XK_adm1.json",
            json2=f"{OSM_BOUNDARIES_DATASET['folder']}/UA_adm1.json",
            json3=f"{OSM_BOUNDARIES_DATASET['folder']}/MD_adm1.json",
            json4=f"{OSM_BOUNDARIES_DATASET['folder']}/BA_adm1.json",
            zip_file=f"{OSM_BOUNDARIES_DATASET['folder']}/osm_boundaries.zip",
        run:
            output_folder = Path(output["zip_file"]).parent
            copy2(input[0], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)


if (
    GEOTHERMAL_HEAT_UTILISATION_POTENTIALS_DATASET := dataset_version(
        "geothermal_heat_utilisation_potentials"
    )
)["source"] in ["primary", "archive"]:

    rule retrieve_geothermal_heat_utilisation_potentials:
        """Retrieves Fraunhofer ISI geothermal heat utilisation potentials."""
        input:
            isi_heat_potentials=storage(
                GEOTHERMAL_HEAT_UTILISATION_POTENTIALS_DATASET["url"]
            ),
        output:
            isi_heat_potentials=f"{GEOTHERMAL_HEAT_UTILISATION_POTENTIALS_DATASET['folder']}/isi_heat_utilisation_potentials.xlsx",
        log:
            "logs/retrieve_geothermal_heat_utilisation_potentials.log",
        retries: 2
        threads: 1
        run:
            copy2(input["isi_heat_potentials"], output["isi_heat_potentials"])


if (LAU_REGIONS_DATASET := dataset_version("lau_regions"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_lau_regions:
        """Retrieves Eurostat Local Administrative Unit region shapes."""
        input:
            lau_regions=storage(LAU_REGIONS_DATASET["url"]),
        output:
            zip=f"{LAU_REGIONS_DATASET['folder']}/lau_regions.zip",
        log:
            "logs/retrieve_lau_regions.log",
        retries: 2
        threads: 1
        run:
            copy2(input["lau_regions"], output["zip"])

    rule retrieve_seawater_temperature:
        """Downloads daily seawater temperature for one year from the Copernicus Marine Service."""
        output:
            seawater_temperature="data/seawater_temperature_{year}.nc",
        log:
            "logs/retrieve_seawater_temperature_{year}.log",
        resources:
            mem_mb=10000,
        params:
            default_cutout=config_provider("atlite", "default_cutout"),
            test_data_url=dataset_version("seawater_temperature")["url"],
        script:
            scripts("retrieve_seawater_temperature.py")

    rule retrieve_hera_data_test_cutout:
        """Retrieves and unpacks a small HERA test dataset for Belgium from Zenodo."""
        input:
            hera_data_url=storage(
                f"https://zenodo.org/records/15828866/files/hera_be_2013-03-01_to_2013-03-08.zip"
            ),
        output:
            river_discharge=f"data/hera_be_2013-03-01_to_2013-03-08/river_discharge_be_2013-03-01_to_2013-03-08.nc",
            ambient_temperature=f"data/hera_be_2013-03-01_to_2013-03-08/ambient_temp_be_2013-03-01_to_2013-03-08.nc",
        log:
            "logs/retrieve_hera_data_test_cutout.log",
        retries: 2
        resources:
            mem_mb=10000,
        params:
            folder="data",
        run:
            unpack_archive(input[0], params.folder)

    rule retrieve_hera_data:
        """Downloads HERA river discharge and air temperature for one year from the JRC."""
        input:
            river_discharge=storage(
                "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-EFAS/HERA/VER1-0/Data/NetCDF/river_discharge/dis.HERA{year}.nc"
            ),
            ambient_temperature=storage(
                "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-EFAS/HERA/VER1-0/Data/NetCDF/climate_inputs/ta6/ta6_{year}.nc"
            ),
        output:
            river_discharge="data/hera_{year}/river_discharge_{year}.nc",
            ambient_temperature="data/hera_{year}/ambient_temp_{year}.nc",
        log:
            "logs/retrieve_hera_data_{year}.log",
        retries: 2
        resources:
            mem_mb=10000,
        params:
            snapshot_year="{year}",
        run:
            move(input.river_discharge, output.river_discharge)
            move(input.ambient_temperature, output.ambient_temperature)


if (JRC_ARDECO_DATASET := dataset_version("jrc_ardeco"))["source"] in [
    "primary",
]:

    rule retrieve_jrc_ardeco:
        """Downloads JRC ARDECO regional GDP and population tables from the ARDECO API."""
        input:
            ardeco_gdp=storage(
                f"{JRC_ARDECO_DATASET['url']}/SUVGDP?versions=2021&unit=EUR&format=csv-table"
            ),
            ardeco_pop=storage(
                f"{JRC_ARDECO_DATASET['url']}/SNPTD?versions=2021&unit=EUR&format=csv-table"
            ),
        output:
            ardeco_gdp=f"{JRC_ARDECO_DATASET['folder']}/ARDECO-SUVGDP.2021.table.csv",
            ardeco_pop=f"{JRC_ARDECO_DATASET['folder']}/ARDECO-SNPTD.2021.table.csv",
        run:
            for key in input.keys():
                copy2(input[key], output[key])

elif (JRC_ARDECO_DATASET := dataset_version("jrc_ardeco"))["source"] in ["archive"]:

    rule retrieve_jrc_ardeco:
        """Retrieves JRC ARDECO regional GDP and population tables from the PyPSA data archive."""
        input:
            ardeco_gdp=storage(
                f"{JRC_ARDECO_DATASET['url']}/ARDECO-SUVGDP.2021.table.csv"
            ),
            ardeco_pop=storage(
                f"{JRC_ARDECO_DATASET['url']}/ARDECO-SNPTD.2021.table.csv"
            ),
        output:
            ardeco_gdp=f"{JRC_ARDECO_DATASET['folder']}/ARDECO-SUVGDP.2021.table.csv",
            ardeco_pop=f"{JRC_ARDECO_DATASET['folder']}/ARDECO-SNPTD.2021.table.csv",
        run:
            for key in input.keys():
                copy2(input[key], output[key])


if (AQUIFER_DATA_DATASET := dataset_version("aquifer_data"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_aquifer_data_bgr:
        """Retrieves and unpacks the BGR International Hydrogeological Map of Europe aquifer shapes."""
        input:
            zip_file=storage(AQUIFER_DATA_DATASET["url"]),
        output:
            zip_file=f"{AQUIFER_DATA_DATASET['folder']}/ihme1500_aquif_ec4060_v12_poly.zip",
            aquifer_shapes=expand(
                f"{AQUIFER_DATA_DATASET['folder']}/IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.{{ext}}",
                ext=[
                    "shp",
                    "shx",
                    "dbf",
                    "cpg",
                    "prj",
                    "sbn",
                    "sbx",
                ],
            ),
        run:
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(
                output["zip_file"],
                AQUIFER_DATA_DATASET["folder"],
            )


if (DH_AREAS_DATASET := dataset_version("dh_areas"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_dh_areas:
        """Retrieves Fraunhofer ISI district heating area shapes."""
        input:
            dh_areas=storage(DH_AREAS_DATASET["url"]),
        output:
            dh_areas=f"{DH_AREAS_DATASET['folder']}/dh_areas.gpkg",
        log:
            "logs/retrieve_dh_areas.log",
        run:
            copy2(input["dh_areas"], output["dh_areas"])


if (MOBILITY_PROFILES_DATASET := dataset_version("mobility_profiles"))["source"] in [
    "archive"
]:

    rule retrieve_mobility_profiles:
        """Retrieves German vehicle activity profiles derived from BASt traffic counts."""
        input:
            kfz=storage(MOBILITY_PROFILES_DATASET["url"] + "/kfz.csv"),
            pkw=storage(MOBILITY_PROFILES_DATASET["url"] + "/pkw.csv"),
        output:
            kfz=f"{MOBILITY_PROFILES_DATASET['folder']}/kfz.csv",
            pkw=f"{MOBILITY_PROFILES_DATASET['folder']}/pkw.csv",
        log:
            "logs/retrieve_mobility_profiles.log",
        benchmark:
            "benchmarks/retrieve_mobility_profiles"
        threads: 1
        resources:
            mem_mb=1000,
        run:
            copy2(input["kfz"], output["kfz"])
            copy2(input["pkw"], output["pkw"])
