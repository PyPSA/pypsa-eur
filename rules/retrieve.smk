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
        message:
            "Retrieving Eurostat balances data"
        input:
            zip_file=storage(EUROSTAT_BALANCES_DATASET["url"]),
        output:
            zip_file=f"{EUROSTAT_BALANCES_DATASET['folder']}/balances.zip",
            directory=directory(EUROSTAT_BALANCES_DATASET["folder"]),
        run:
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output["directory"])


if (
    EUROSTAT_HOUSEHOLD_BALANCES_DATASET := dataset_version(
        "eurostat_household_balances"
    )
)["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eurostat_household_balances:
        message:
            "Retrieving Eurostat household balances data"
        input:
            csv=storage(EUROSTAT_HOUSEHOLD_BALANCES_DATASET["url"]),
        output:
            csv=f"{EUROSTAT_HOUSEHOLD_BALANCES_DATASET['folder']}/nrg_d_hhq.csv",
        run:
            copy2(input["csv"], output["csv"])


if (NUTS3_POPULATION_DATASET := dataset_version("nuts3_population"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_nuts3_population:
        message:
            "Retrieving NUTS3 population data"
        input:
            gz=storage(NUTS3_POPULATION_DATASET["url"]),
        output:
            gz=f"{NUTS3_POPULATION_DATASET['folder']}/nama_10r_3popgdp.tsv.gz",
        retries: 2
        run:
            copy2(input["gz"], output["gz"])


if (CORINE_DATASET := dataset_version("corine"))["source"] in ["archive"]:

    rule retrieve_corine:
        message:
            "Retrieving Corine land cover data"
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
        message:
            "Retrieving Corine land cover data"
        params:
            apikey=os.environ.get("CORINE_API_TOKEN", ""),
        output:
            zip=f"{CORINE_DATASET['folder']}/corine.zip",
            tif_file=f"{CORINE_DATASET['folder']}/corine.tif",
        log:
            logs("retrieve_corine_primary.log"),
        resources:
            mem_mb=1000,
        retries: 2
        script:
            "../scripts/retrieve_corine_dataset_primary.py"


if (H2_SALT_CAVERNS_DATASET := dataset_version("h2_salt_caverns"))["source"] in [
    "archive"
]:

    rule retrieve_h2_salt_caverns:
        message:
            "Retrieving H2 salt caverns data"
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
        message:
            "Retrieving GDP per capita data"
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
        message:
            "Retrieving population count data"
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
        message:
            "Retrieving GHG emissions data"
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
        message:
            "Retrieving GEBCO bathymetry data"
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
        message:
            "Retrieving attributed ports data"
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
        message:
            "Retrieving JRC IDEES data"
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
        message:
            "Retrieving EU NUTS 2013 data"
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
        message:
            "Retrieving EU NUTS 2021 data"
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
        input:
            geojson=storage(BIDDING_ZONES_ELECTRICITYMAPS_DATASET["url"]),
        output:
            geojson=f"{BIDDING_ZONES_ELECTRICITYMAPS_DATASET['folder']}/bidding_zones_electricitymaps.geojson",
        log:
            "logs/retrieve_bidding_zones_electricitymaps.log",
        resources:
            mem_mb=1000,
        retries: 2
        run:
            copy2(input["geojson"], output["geojson"])


if (BIDDING_ZONES_ENTSOEPY_DATASET := dataset_version("bidding_zones_entsoepy"))[
    "source"
] in ["primary", "archive"]:

    rule retrieve_bidding_zones_entsoepy:
        output:
            geojson=f"{BIDDING_ZONES_ENTSOEPY_DATASET['folder']}/bidding_zones_entsoepy.geojson",
        log:
            "logs/retrieve_bidding_zones_entsoepy.log",
        resources:
            mem_mb=1000,
        retries: 2
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
        message:
            "Retrieving cutout data for {wildcards.cutout}"
        input:
            storage(CUTOUT_DATASET["url"] + "/{cutout}.nc"),
        output:
            CUTOUT_DATASET["folder"] + "/{cutout}.nc",
        log:
            "logs/retrieve_cutout/{cutout}.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            copy2(input[0], output[0])


if (COUNTRY_RUNOFF_DATASET := dataset_version("country_runoff"))["source"] in [
    "archive"
]:

    rule retrieve_country_runoff:
        message:
            "Retrieving country runoff data"
        input:
            storage(COUNTRY_RUNOFF_DATASET["url"]),
        output:
            era5_runoff=f"{COUNTRY_RUNOFF_DATASET['folder']}/era5-runoff-per-country.csv",
        run:
            copy2(input[0], output[0])


if (COUNTRY_HDD_DATASET := dataset_version("country_hdd"))["source"] in ["archive"]:

    rule retrieve_country_hdd:
        message:
            "Retrieving country heating degree days data"
        input:
            storage(COUNTRY_HDD_DATASET["url"]),
        output:
            era5_runoff=f"{COUNTRY_HDD_DATASET['folder']}/era5-HDD-per-country.csv",
        run:
            copy2(input[0], output[0])


if (COSTS_DATASET := dataset_version("costs"))["source"] in [
    "primary",
]:

    rule retrieve_cost_data:
        message:
            "Retrieving cost data for {wildcards.planning_horizons}"
        input:
            costs=storage(COSTS_DATASET["url"] + "/costs_{planning_horizons}.csv"),
        output:
            costs=COSTS_DATASET["folder"] + "/costs_{planning_horizons}.csv",
        run:
            copy2(input["costs"], output["costs"])


if (POWERPLANTS_DATASET := dataset_version("powerplants"))["source"] in [
    "primary",
]:

    rule retrieve_powerplants:
        message:
            "Retrieving powerplants data"
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
        message:
            "Retrieving SciGRID gas infrastructure data"
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


rule retrieve_electricity_demand:
    message:
        "Retrieving electricity demand data"
    params:
        versions=["2019-06-05", "2020-10-06"],
    output:
        "data/electricity_demand_raw.csv",
    log:
        "logs/retrieve_electricity_demand.log",
    resources:
        mem_mb=5000,
    retries: 2
    script:
        "../scripts/retrieve_electricity_demand.py"


if (
    SYNTHETIC_ELECTRICITY_DEMAND_DATASET := dataset_version(
        "synthetic_electricity_demand"
    )
)["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_synthetic_electricity_demand:
        message:
            "Retrieving synthetic electricity demand data"
        input:
            csv=storage(SYNTHETIC_ELECTRICITY_DEMAND_DATASET["url"]),
        output:
            csv=f"{SYNTHETIC_ELECTRICITY_DEMAND_DATASET['folder']}/load_synthetic_raw.csv",
        retries: 2
        run:
            copy2(input["csv"], output["csv"])


if (SHIP_RASTER_DATASET := dataset_version("ship_raster"))["source"] in [
    "archive",
    "primary",
]:

    rule retrieve_ship_raster:
        message:
            "Retrieving shipping raster data"
        input:
            zip_file=storage(SHIP_RASTER_DATASET["url"]),
        output:
            zip_file=f"{SHIP_RASTER_DATASET['folder']}/shipdensity_global.zip",
        log:
            "logs/retrieve_ship_raster.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            copy2(input["zip_file"], output["zip_file"])


if (ENSPRESO_BIOMASS_DATASET := dataset_version("enspreso_biomass"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_enspreso_biomass:
        message:
            "Retrieving ENSPRESO biomass data"
        input:
            xlsx=storage(ENSPRESO_BIOMASS_DATASET["url"]),
        output:
            xlsx=f"{ENSPRESO_BIOMASS_DATASET['folder']}/ENSPRESO_BIOMASS.xlsx",
        retries: 1
        run:
            copy2(input["xlsx"], output["xlsx"])


if (HOTMAPS_INDUSTRIAL_SITES := dataset_version("hotmaps_industrial_sites"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_hotmaps_industrial_sites:
        message:
            "Retrieving Hotmaps industrial sites"
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
        message:
            "Retrieving nitrogen statistics data"
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
        message:
            "Retrieving Copernicus land cover data"
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
        message:
            "Retrieving LUISA land cover data"
        input:
            tif=storage(LUISA_LAND_COVER_DATASET["url"]),
        output:
            tif=f"{LUISA_LAND_COVER_DATASET['folder']}/LUISA_basemap_020321_50m.tif",
        run:
            copy2(input["tif"], output["tif"])


if (EEZ_DATASET := dataset_version("eez"))["source"] in ["primary"]:

    rule retrieve_eez:
        message:
            "Retrieving EEZ data"
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
        message:
            "Retrieving EEZ data"
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
        message:
            "Retrieving World Bank urban population data"
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
        message:
            "Retrieving CO2STOP data"
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
        message:
            "Retrieving GEM Europe Gas Tracker data"
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
        message:
            "Retrieving GEM Global Steel Plant Tracker data"
        input:
            xlsx=storage(GEM_GSPT_DATASET["url"]),
        output:
            xlsx=f"{GEM_GSPT_DATASET['folder']}/Global-Steel-Plant-Tracker.xlsx",
        run:
            copy2(input["xlsx"], output["xlsx"])


if (BFS_ROAD_VEHICLE_STOCK_DATASET := dataset_version("bfs_road_vehicle_stock"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_bfs_road_vehicle_stock:
        message:
            "Retrieving BFS road vehicle stock data"
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
        message:
            "Retrieving BFS GDP and population data"
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
        message:
            "Downloading protected area database from WDPA"
        input:
            zip_file=storage(get_wdpa_url(WDPA_DATASET)),
        output:
            zip_file=f"{WDPA_DATASET['folder']}/WDPA_shp.zip",
            gpkg=f"{WDPA_DATASET['folder']}/WDPA.gpkg",
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
        message:
            "Downloading Marine protected area database from WDPA"
        input:
            zip_file=storage(get_wdpa_url(WDPA_MARINE_DATASET)),
        output:
            zip_file=f"{WDPA_MARINE_DATASET['folder']}/WDPA_WDOECM_marine.zip",
            gpkg=f"{WDPA_MARINE_DATASET['folder']}/WDPA_WDOECM_marine.gpkg",
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



# Versioning not implemented as the dataset is used only for validation
# License - (c) EEX AG, all rights reserved. Personal copy for non-commercial use permitted
rule retrieve_monthly_co2_prices:
    message:
        "Retrieving monthly CO2 prices data for validation"
    input:
        storage(
            "https://public.eex-group.com/eex/eua-auction-report/emission-spot-primary-market-auction-report-2019-data.xls",
        ),
    output:
        "data/validation/emission-spot-primary-market-auction-report-2019-data.xls",
    log:
        "logs/retrieve_monthly_co2_prices.log",
    resources:
        mem_mb=5000,
    retries: 2
    run:
        copy2(input[0], output[0])


# Versioning not implemented as the dataset is used only for validation
# License - custom; no restrictions on use and redistribution, attribution required
rule retrieve_monthly_fuel_prices:
    message:
        "Retrieving monthly fuel prices data for validation"
    output:
        "data/validation/energy-price-trends-xlsx-5619002.xlsx",
    log:
        "logs/retrieve_monthly_fuel_prices.log",
    resources:
        mem_mb=5000,
    retries: 2
    script:
        "../scripts/retrieve_monthly_fuel_prices.py"


if (TYDNP_DATASET := dataset_version("tyndp"))["source"] in ["primary", "archive"]:

    rule retrieve_tyndp:
        message:
            "Retrieving TYNDP network topology data"
        input:
            line_data=storage(TYDNP_DATASET["url"] + "/Line-data.zip"),
            nodes=storage(TYDNP_DATASET["url"] + "/Nodes.zip"),
        output:
            line_data_zip=f"{TYDNP_DATASET['folder']}/Line-data.zip",
            nodes_zip=f"{TYDNP_DATASET['folder']}/Nodes.zip",
            reference_grid=f"{TYDNP_DATASET['folder']}/Line data/ReferenceGrid_Electricity.xlsx",
            nodes=f"{TYDNP_DATASET['folder']}/Nodes/LIST OF NODES.xlsx",
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



if OSM_DATASET["source"] in ["archive"]:

    OSM_ARCHIVE_FILES = [
        "buses.csv",
        "converters.csv",
        "lines.csv",
        "links.csv",
        "transformers.csv",
        # Newer versions include the additional map.html file for visualisation
        *(["map.html"] if float(OSM_DATASET["version"]) >= 0.6 else []),
    ]

    rule retrieve_osm_archive:
        message:
            "Retrieving OSM archive data"
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


elif OSM_DATASET["source"] == "build":

    OSM_RAW_JSON = [
        "cables_way.json",
        "lines_way.json",
        "routes_relation.json",
        "substations_way.json",
        "substations_relation.json",
    ]

    rule retrieve_osm_data_raw:
        message:
            "Retrieving OSM raw data for {wildcards.country}"
        params:
            overpass_api=config_provider("overpass_api"),
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
        script:
            "../scripts/retrieve_osm_data.py"

    rule retrieve_osm_data_raw_all:
        input:
            expand(
                f"{OSM_DATASET['folder']}/{{country}}/{{file}}",
                country=config_provider("countries"),
                file=OSM_RAW_JSON,
            ),


if (NATURA_DATASET := dataset_version("natura"))["source"] in ["archive"]:

    rule retrieve_natura:
        message:
            "Retrieving Natura 2000 raster data"
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
        message:
            "Building Natura 2000 raster data"
        input:
            online=storage(NATURA_DATASET["url"]),
            cutout=lambda w: input_cutout(w),
        output:
            zip=f"{NATURA_DATASET['folder']}/raw/natura.zip",
            raw=directory(f"{NATURA_DATASET['folder']}/raw"),
            raster=f"{NATURA_DATASET['folder']}/natura.tiff",
        resources:
            mem_mb=5000,
        log:
            "logs/build_natura.log",
        script:
            "../scripts/build_natura.py"


if (OSM_BOUNDARIES_DATASET := dataset_version("osm_boundaries"))["source"] in [
    "primary"
]:

    rule retrieve_osm_boundaries:
        message:
            "Retrieving OSM admin boundaries for {wildcards.country}"
        output:
            json=f"{OSM_BOUNDARIES_DATASET['folder']}/{country}_adm1.json",
        log:
            "logs/retrieve_osm_boundaries_{country}_adm1.log",
        threads: 1
        script:
            "../scripts/retrieve_osm_boundaries.py"

elif (OSM_BOUNDARIES_DATASET := dataset_version("osm_boundaries"))["source"] in [
    "archive"
]:

    rule retrieve_osm_boundaries:
        message:
            "Retrieving OSM admin boundaries data"
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
        message:
            "Retrieving geothermal heat utilisation potentials"
        input:
            isi_heat_potentials=storage(
                GEOTHERMAL_HEAT_UTILISATION_POTENTIALS_DATASET["url"]
            ),
        output:
            isi_heat_potentials=f"{GEOTHERMAL_HEAT_UTILISATION_POTENTIALS_DATASET['folder']}/isi_heat_utilisation_potentials.xlsx",
        log:
            "logs/retrieve_geothermal_heat_utilisation_potentials.log",
        threads: 1
        retries: 2
        run:
            copy2(input["isi_heat_potentials"], output["isi_heat_potentials"])


if (LAU_REGIONS_DATASET := dataset_version("lau_regions"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_lau_regions:
        message:
            "Retrieving Local Administrative Units and Administation Unit regions"
        input:
            lau_regions=storage(LAU_REGIONS_DATASET["url"]),
        output:
            zip=f"{LAU_REGIONS_DATASET['folder']}/lau_regions.zip",
        log:
            "logs/retrieve_lau_regions.log",
        threads: 1
        retries: 2
        run:
            copy2(input["lau_regions"], output["zip"])

    rule retrieve_seawater_temperature:
        message:
            "Retrieving seawater temperature data for {wildcards.year}"
        params:
            default_cutout=config_provider("atlite", "default_cutout"),
            test_data_url=dataset_version("seawater_temperature")["url"],
        output:
            seawater_temperature="data/seawater_temperature_{year}.nc",
        log:
            "logs/retrieve_seawater_temperature_{year}.log",
        resources:
            mem_mb=10000,
        script:
            "../scripts/retrieve_seawater_temperature.py"

    rule retrieve_hera_data_test_cutout:
        message:
            "Retrieving HERA test cutout data"
        input:
            hera_data_url=storage(
                f"https://zenodo.org/records/15828866/files/hera_be_2013-03-01_to_2013-03-08.zip"
            ),
        output:
            river_discharge=f"data/hera_be_2013-03-01_to_2013-03-08/river_discharge_be_2013-03-01_to_2013-03-08.nc",
            ambient_temperature=f"data/hera_be_2013-03-01_to_2013-03-08/ambient_temp_be_2013-03-01_to_2013-03-08.nc",
        params:
            folder="data",
        log:
            "logs/retrieve_hera_data_test_cutout.log",
        resources:
            mem_mb=10000,
        retries: 2
        run:
            unpack_archive(input[0], params.folder)

    rule retrieve_hera_data:
        message:
            "Retrieving HERA data for {wildcards.year}"
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
        params:
            snapshot_year="{year}",
        log:
            "logs/retrieve_hera_data_{year}.log",
        resources:
            mem_mb=10000,
        retries: 2
        run:
            move(input.river_discharge, output.river_discharge)
            move(input.ambient_temperature, output.ambient_temperature)


if (JRC_ARDECO_DATASET := dataset_version("jrc_ardeco"))["source"] in [
    "primary",
]:

    rule retrieve_jrc_ardeco:
        message:
            "Retrieving JRC ARDECO data"
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
        message:
            "Retrieving JRC ARDECO data"
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
        message:
            "Retrieving BGR aquifer data"
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
        message:
            "Retrieving District Heating areas"
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
        message:
            "Retrieving mobility profiles data"
        input:
            kfz=storage(MOBILITY_PROFILES_DATASET["url"] + "/kfz.csv"),
            pkw=storage(MOBILITY_PROFILES_DATASET["url"] + "/pkw.csv"),
        output:
            kfz=f"{MOBILITY_PROFILES_DATASET['folder']}/kfz.csv",
            pkw=f"{MOBILITY_PROFILES_DATASET['folder']}/pkw.csv",
        threads: 1
        resources:
            mem_mb=1000,
        log:
            "logs/retrieve_mobility_profiles.log",
        benchmark:
            "benchmarks/retrieve_mobility_profiles"
        run:
            copy2(input["kfz"], output["kfz"])
            copy2(input["pkw"], output["pkw"])
