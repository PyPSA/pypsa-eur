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
storage:
    provider="http",
    keep_local=True,
    retrieve=True,
    retries=2,


if config["enable"].get("retrieve", "auto") == "auto":
    config["enable"]["retrieve"] = has_internet_access()

if config["enable"]["retrieve"] is False:
    print("Datafile downloads disabled in config[retrieve] or no internet access.")


if (EUROSTAT_BALANCES_DATASET := dataset_version("eurostat_balances"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eurostat_balances:
        input:
            zip_file=storage(
                EUROSTAT_BALANCES_DATASET["url"],
            ),
        output:
            zip_file=f"{EUROSTAT_BALANCES_DATASET['folder']}/balances.zip",
            directory=directory(f"{EUROSTAT_BALANCES_DATASET['folder']}"),
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
        input:
            csv=storage(
                EUROSTAT_HOUSEHOLD_BALANCES_DATASET["url"],
            ),
        output:
            csv=f"{EUROSTAT_HOUSEHOLD_BALANCES_DATASET['folder']}/nrg_d_hhq.csv",
        run:
            copy2(input["csv"], output["csv"])

if (
    NUTS3_POPULATION_DATASET := dataset_version("nuts3_population")
)["source"] in [
    "primary", 
    "archive"
]:

    rule retrieve_nuts3_population:
        input:
            storage(
                NUTS3_POPULATION_DATASET["url"],
            ),
        output:
            gz=f"{NUTS3_POPULATION_DATASET["folder"]}/nama_10r_3popgdp.tsv.gz",
        retries: 2
        run:
            move(input[0], output[0])



if (CORINE_DATASET := dataset_version("corine"))["source"] in [
    "archive"
]:
    rule retrieve_corine:
        params:
            url = CORINE_DATASET["url"],
        output:
            zip=f"{CORINE_DATASET["folder"]}/corine.zip",
            directory=directory(f"{CORINE_DATASET["folder"]}"),
            tif_file=f"{CORINE_DATASET["folder"]}/corine/g250_clc06_V18_5.tif",
        run:
            response = requests.get(params["url"])
            with open(output.zip, "wb") as f:
                f.write(response.content)

            output_folder = Path(output["zip"]).parent
            unpack_archive(output.zip, output_folder)

elif (CORINE_DATASET := dataset_version("corine"))["source"] in [
    "primary"
]:
    rule retrieve_corine:
        params:
            apikey = config["secrets"]["corine"],
        output:
            zip=f"{CORINE_DATASET["folder"]}/corine.zip",
            tif_file_path=f"{CORINE_DATASET["folder"]}/Results/u2018_clc2012_v2020_20u1_raster100m/DATA/U2018_CLC2012_V2020_20u1.tif",
        log:
            logs("retrieve_corine_primary.log"),
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_corine_dataset_primary.py"


if (H2_SALT_CAVERNS_DATASET := dataset_version("h2_salt_caverns"))["source"] in [
    "archive"
]:
    rule retrieve_h2_salt_caverns:
        input: 
            storage(
                H2_SALT_CAVERNS_DATASET["url"],
            ),
        output:
            geojson=f"{H2_SALT_CAVERNS_DATASET["folder"]}/h2_salt_caverns_GWh_per_sqkm.geojson",
        retries: 2,
        run:
            move(input[0], output[0])


if (
    GDP_PER_CAPITA_DATASET := dataset_version("gdp_per_capita")
)["source"] in ["archive"]:
    rule retrieve_gdp_per_capita:
        input:
            storage(
                GDP_PER_CAPITA_DATASET["url"],
            ),
        output:
            gdp = f"{GDP_PER_CAPITA_DATASET["folder"]}/GDP_per_capita_PPP_1990_2015_v2.nc",
        retries: 2
        run:
            move(input[0], output[0])


if (
    POPULATION_COUNT_DATASET := dataset_version("population_count")
)["source"] in ["archive", "primary"]:

    rule retrieve_population_count:
        input:
            storage(
                POPULATION_COUNT_DATASET["url"],
            ),
        output:
            tif = f"{POPULATION_COUNT_DATASET["folder"]}/ppp_2019_1km_Aggregated.tif",
        retries: 2
        run:
            move(input[0], output[0])

            if POPULATION_COUNT_DATASET["source"] == 'primary':
                import xarray as xr
                import rioxarray as rio
                file_path = output[0]
                ds = xr.open_dataarray(file_path)
                ds_reqd = ds.sel(x=slice(15.55, 40.41), y=slice(52.49,41.72))
                ds_reqd.rio.to_raster(file_path)
            

if (
    GHG_EMISSIONS_DATASET := dataset_version("ghg_emissions")
)["source"] in ["archive", "primary"]:

    rule retrieve_ghg_emissions:
        params:
            url = GHG_EMISSIONS_DATASET["url"],
        input:
            storage(
                GHG_EMISSIONS_DATASET["url"],
            ),
        output:
            csv = f"{GHG_EMISSIONS_DATASET["folder"]}/UNFCCC_v23.csv",
            zip=f"{GHG_EMISSIONS_DATASET["folder"]}/UNFCCC_v23.csv.zip" if GHG_EMISSIONS_DATASET["source"] == "primary" else [],
            directory=directory(f"{GHG_EMISSIONS_DATASET["folder"]}") if GHG_EMISSIONS_DATASET["source"] == "primary" else [],
        retries: 2
        run:
            move(input[0], output[0])
            if GHG_EMISSIONS_DATASET["source"] == "primary":
                handle_data_requests(params, output)


if (GEBCO_DATASET := dataset_version("gebco"))["source"] in [
    "archive",
    "primary"
]:
    rule retrieve_gebco:
        params:
            url = GEBCO_DATASET["url"],
        input:
            storage(
                GEBCO_DATASET["url"],
            ) if GEBCO_DATASET["source"] == 'archive' else [],
        output:
            gebco=f"{GEBCO_DATASET["folder"]}/GEBCO_2014_2D.nc",
            zip=f"{GEBCO_DATASET["folder"]}/GEBCO_2014.zip" if GEBCO_DATASET["source"] == 'primary' else [],
            directory=directory(f"{GEBCO_DATASET["folder"]}") if GEBCO_DATASET["source"] == 'primary' else [],
        retries:2,
        run:
            if GEBCO_DATASET["source"] == 'primary':
                import xarray as xr

                response = requests.get(params["url"])
                with open(output.zip, "wb") as f:
                    f.write(response.content)

                output_folder = Path(output["zip"]).parent
                unpack_archive(output.zip, output_folder)

                file_path = output["gebco"]
                ds = xr.open_dataset(file_path)
                ds_reqd = ds.sel(lat=slice(32,73),lon=slice(-21,45))
                ds_reqd.to_netcdf(file_path)
            else:
                move(input[0], output[0])
                validate_checksum(output[0], input[0])


if (ATTRIBUTED_PORTS_DATASET := dataset_version("attributed_ports"))["source"] in [
    "archive",
    "primary"
]:
    rule retrieve_attributed_ports:
        input:
            storage(
                ATTRIBUTED_PORTS_DATASET["url"],
            )
        output:
            json=f"{ATTRIBUTED_PORTS_DATASET["folder"]}/attributed_ports.json",
        retries: 2,
        run:
            move(input[0], output[0])

if (JRC_IDEES_DATASET := dataset_version("jrc_idees"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_jrc_idees:
        params:
            url=JRC_IDEES_DATASET["url"],
        output:
            zip=f"{JRC_IDEES_DATASET["folder"]}/jrc_idees.zip",
            directory=directory(f"{JRC_IDEES_DATASET["folder"]}"),
        run:
            response = requests.get(params["url"])
            with open(output.zip, "wb") as f:
                f.write(response.content)

            output_folder = Path(output["zip"]).parent
            unpack_archive(output.zip, output_folder)


if (EU_NUTS2013_DATASET := dataset_version("eu_nuts2013"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eu_nuts_2013:
        input:
            shapes=storage(EU_NUTS2013_DATASET["url"]),
        output:
            zip_file=f"{EU_NUTS2013_DATASET['folder']}/ref-nuts-2013-03m.geojson.zip",
            folder=directory(
                f"{EU_NUTS2013_DATASET["folder"]}/ref-nuts-2013-03m.geojson"
            ),
            shapes_level_3=f"{EU_NUTS2013_DATASET["folder"]}/ref-nuts-2013-03m.geojson/NUTS_RG_03M_2013_4326_LEVL_3.geojson",
            shapes_level_2=f"{EU_NUTS2013_DATASET["folder"]}/ref-nuts-2013-03m.geojson/NUTS_RG_03M_2013_4326_LEVL_2.geojson",
        run:
            copy2(input["shapes"], output["zip_file"])
            unpack_archive(output["zip_file"], Path(output.shapes_level_3).parent)


if (EU_NUTS2021_DATASET := dataset_version("eu_nuts2021"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eu_nuts_2021:
        input:
            shapes=storage(
                EU_NUTS2021_DATASET["url"],
            ),
        output:
            zip_file=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson.zip",
            folder=directory(
                f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson"
            ),
            shapes_level_3=f"{EU_NUTS2021_DATASET["folder"]}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_3.geojson",
            shapes_level_2=f"{EU_NUTS2021_DATASET["folder"]}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_2.geojson",
            shapes_level_1=f"{EU_NUTS2021_DATASET["folder"]}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_1.geojson",
            shapes_level_0=f"{EU_NUTS2021_DATASET["folder"]}/ref-nuts-2021-01m.geojson/NUTS_RG_01M_2021_4326_LEVL_0.geojson",
        run:
            copy2(input["shapes"], output["zip_file"])
            unpack_archive(output["zip_file"], Path(output.shapes_level_3).parent)


if config["enable"]["retrieve"]:

    rule retrieve_bidding_zones:
        output:
            file_entsoepy="data/busshapes/bidding_zones_entsoepy.geojson",
            file_electricitymaps="data/busshapes/bidding_zones_electricitymaps.geojson",
        log:
            "logs/retrieve_bidding_zones.log",
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_bidding_zones.py"


if (CUTOUT_DATASET := dataset_version("cutout"))["source"] in [
    "archive",
]:

    rule retrieve_cutout:
        input:
            storage(
                CUTOUT_DATASET["url"] + "/files/{cutout}.nc",
            ),
        output:
            CUTOUT_DATASET["folder"] / "{cutout}.nc",
        log:
            "logs/retrieve_cutout/{cutout}.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])
            validate_checksum(output[0], input[0])


if (COUNTRY_RUNOFF_DATASET := dataset_version("country_runoff"))["source"] in [
    "archive"
]:

    rule retrieve_country_runoff:
        input:
            storage(COUNTRY_RUNOFF_DATASET["url"]),
        output:
            era5_runoff=COUNTRY_RUNOFF_DATASET["folder"] / "era5-runoff-per-country.csv",
        run:
            move(input[0], output[0])


if (COUNTRY_HDD_DATASET := dataset_version("country_hdd"))["source"] in ["archive"]:

    rule retrieve_country_hdd:
        input:
            storage(COUNTRY_HDD_DATASET["url"]),
        output:
            era5_runoff=COUNTRY_HDD_DATASET["folder"] / "era5-HDD-per-country.csv",
        run:
            move(input[0], output[0])


if (COSTS_DATASET := dataset_version("costs"))["source"] in [
    "primary",
]:

    rule retrieve_cost_data:
        input:
            storage(COSTS_DATASET["url"] + "/costs_{year}.csv"),
        output:
            costs=COSTS_DATASET["folder"] / "costs_{year}.csv",
        run:
            move(input[0], output[0])


if (POWERPLANTS_DATASET := dataset_version("powerplants"))["source"] in [
    "primary",
]:

    rule retrieve_powerplants:
        input:
            storage(POWERPLANTS_DATASET["url"]),
        output:
            powerplants=f"{POWERPLANTS_DATASET['folder']}/powerplants.csv",
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"] and (
    SCIGRID_GAS_DATASET := dataset_version("scigrid_gas")
)["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_gas_infrastructure_data:
        params:
            url=SCIGRID_GAS_DATASET['url'],
        output:
            zip=f"{SCIGRID_GAS_DATASET["folder"]}/jrc_idees.zip",
            directory=directory(f"{SCIGRID_GAS_DATASET["folder"]}"),
        run:

            response = requests.get(params["url"])
            with open(output.zip, "wb") as f:
                f.write(response.content)

            output_folder = Path(output["zip"]).parent
            unpack_archive(output.zip, output_folder)



if config["enable"]["retrieve"]:

    rule retrieve_electricity_demand:
        params:
            versions=["2019-06-05", "2020-10-06"],
        output:
            "data/electricity_demand_raw.csv",
        log:
            "logs/retrieve_electricity_demand.log",
        resources:
            mem_mb=5000,
        retries: 2
        conda:
            "../envs/environment.yaml"
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
        input:
            storage(
                SYNTHETIC_ELECTRICITY_DEMAND_DATASET["url"],
                keep_local=True,
            ),
        output:
            csv=f"{SYNTHETIC_ELECTRICITY_DEMAND_DATASET["folder"]}/load_synthetic_raw.csv",
        retries: 2
        run:
            move(input[0], output[0])


if (
    SHIP_RASTER_DATASET := dataset_version("ship_raster")
)["source"] in ["archive","primary"]:

    rule retrieve_ship_raster:
        input:
            storage(
                SHIP_RASTER_DATASET["url"],
                keep_local=True,
            ),
        output:
            zip_file=f"{SHIP_RASTER_DATASET["folder"]}/shipdensity_global.zip",
        log:
            "logs/retrieve_ship_raster.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


if (ENSPRESO_BIOMASS_DATASET := dataset_version("enspreso_biomass"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_enspreso_biomass:
        input:
            xlsx=storage(
                ENSPRESO_BIOMASS_DATASET["url"],
            ),
        output:
            xlsx=f"{ENSPRESO_BIOMASS_DATASET["folder"]}/ENSPRESO_BIOMASS.xlsx",
        retries: 1
        run:
            move(input["xlsx"], output["xlsx"])


if (HOTMAPS_INDUSTRIAL_SITES := dataset_version("hotmaps_industrial_sites"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_hotmaps_industrial_sites:
        input:
            storage(
                HOTMAPS_INDUSTRIAL_SITES["url"],
            ),
        output:
            csv=f"{HOTMAPS_INDUSTRIAL_SITES["folder"]}/Industrial_Database.csv",
        retries: 1
        run:
            move(input[0], output[0])


if (NITROGEN_STATISTICS_DATASET := dataset_version("nitrogen_statistics"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_nitrogen_statistics:
        input:
            storage(
                NITROGEN_STATISTICS_DATASET["url"],
            ),
        output:
            xlsx=f"{NITROGEN_STATISTICS_DATASET['folder']}/nitro-ert.xlsx",
        retries: 1
        run:
            move(input[0], output[0])


if (
    COPERNICUS_LAND_COVER_DATASET := dataset_version("copernicus_land_cover")
)["source"] in ["primary", "archive"]:

    # Downloading Copernicus Global Land Cover for land cover and land use:
    # Website: https://land.copernicus.eu/global/products/lc
    rule download_copernicus_land_cover:
        input:
            storage(
                COPERNICUS_LAND_COVER_DATASET["url"],
            ),
        output:
            tif = f"{COPERNICUS_LAND_COVER_DATASET["folder"]}/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        run:
            move(input[0], output[0])
            validate_checksum(output[0], input[0])


if (LUISA_LAND_COVER_DATASET := dataset_version("luisa_land_cover"))["source"] in [
    "primary",
    "archive",
]:

    # Downloading LUISA Base Map for land cover and land use:
    # Website: https://ec.europa.eu/jrc/en/luisa
    rule retrieve_luisa_land_cover:
        input:
            storage(
                LUISA_LAND_COVER_DATASET["url"],
            ),
        output:
            tif=f"{LUISA_LAND_COVER_DATASET["folder"]}/LUISA_basemap_020321_50m.tif",
        run:
            move(input[0], output[0])


if (EEZ_DATASET := dataset_version("eez"))[
    "source"
] in ["primary", "archive"]:

    rule retrieve_eez:
        params:
            zip=f"{EEZ_DATASET["folder"]}/World_EEZ_{EEZ_DATASET["version"]}_LR.zip",
        output:
            gpkg=f"{EEZ_DATASET["folder"]}/World_EEZ_{EEZ_DATASET["version"]}_LR/eez_{EEZ_DATASET["version"].split("_")[0]}_lowres.gpkg",
        run:
            from uuid import uuid4

            if EEZ_DATASET["source"] == "primary":

                name = str(uuid4())[:8]
                org = str(uuid4())[:8]

                response = requests.post(
                    f"{EEZ_DATASET["url"]}",
                    params={"name": f"World_EEZ_{EEZ_DATASET["version"]}_LR.zip"},
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

            elif EEZ_DATASET["source"] == "archive":
                response = requests.get(f"{EEZ_DATASET["url"]}")

            with open(params["zip_file"], "wb") as f:
                f.write(response.content)
            output_folder = Path(output.gpkg).parent.parent
            unpack_archive(params["zip_file"], output_folder)
            os.remove(params["zip_file"])



if (WB_URB_POP_DATASET := dataset_version("worldbank_urban_population"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_worldbank_urban_population:
        params:
            url=WB_URB_POP_DATASET["url"],
        output:
            zip=f"{WB_URB_POP_DATASET['folder']}/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.zip",
            csv=f"{WB_URB_POP_DATASET['folder']}/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.csv",
        run:
            response = requests.get(
                params["url"],
            )

            with open(output["zip"], "wb") as f:
                f.write(response.content)
            output_folder = Path(output["zip"]).parent
            unpack_archive(output["zip"], output_folder)

            for f in os.listdir(output_folder):
                if f.startswith(
                    "API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2_"
                ) and f.endswith(".csv"):
                    os.rename(os.path.join(output_folder, f), output.csv)
                    break
            os.remove(params["zip_file"])



if (CO2STOP_DATASET := dataset_version("co2stop"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_co2stop:
        input:
            zip_file=storage(
                CO2STOP_DATASET["url"],
            ),
        output:
            zip_file=f"{CO2STOP_DATASET['folder']}/co2jrc_openformats.zip",
            storage_table=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Storage_Units.csv",
            storage_map=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_Polygons Data/StorageUnits_March13.kml",
            traps_table1=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps.csv",
            traps_table2=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps_Temp.csv",
            traps_table3=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps1.csv",
            traps_map=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_Polygons Data/DaughterUnits_March13.kml",
        run:
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], CO2STOP_DATASET["folder"])


if (GEM_EUROPE_GAS_TRACKER_DATASET := dataset_version("gem_europe_gas_tracker"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_gem_europe_gas_tracker:
        input:
            xlsx=storage(GEM_EUROPE_GAS_TRACKER_DATASET["url"]),
        output:
            xlsx="data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        run:
            move(input["xlsx"], output["xlsx"])


if (GEM_GSPT_DATASET := dataset_version("gem_gspt"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_gem_steel_plant_tracker:
        input:
            xlsx=storage(GEM_GSPT_DATASET["url"]),
        output:
            xlsx=f"{GEM_GSPT_DATASET['folder']}/Global-Steel-Plant-Tracker.xlsx",
        run:
            move(input["xlsx"], output["xlsx"])


if (BFS_ROAD_VEHICLE_STOCK_DATASET := dataset_version("bfs_road_vehicle_stock"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_bfs_road_vehicle_stock:
        input:
            csv=storage(BFS_ROAD_VEHICLE_STOCK_DATASET["url"]),
        output:
            csv=f"{BFS_ROAD_VEHICLE_STOCK_DATASET['folder']}/vehicle_stock.csv",
        run:
            move(input["csv"], output["csv"])


if (BFS_GDP_AND_POPULATION_DATASET := dataset_version("bfs_gdp_and_population"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_bfs_gdp_and_population:
        input:
            xlsx=storage(BFS_GDP_AND_POPULATION_DATASET["url"]),
        output:
            xlsx=f"{BFS_GDP_AND_POPULATION_DATASET['folder']}/gdp_and_population.xlsx",
        run:
            move(input["xlsx"], output["xlsx"])


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
        input:
            zip_file=storage(get_wdpa_url(WDPA_DATASET)),
        output:
            zip_file=f"{WDPA_DATASET['folder']}/WDPA_shp.zip",
            gpkg=f"{WDPA_DATASET['folder']}/WDPA.gpkg",
        run:
            output_folder = Path(output["zip_file"]).parent
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)

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
        input:
            zip_file=storage(
                get_wdpa_url(WDPA_MARINE_DATASET),
            ),
        output:
            zip_file=f"{WDPA_MARINE_DATASET['folder']}/WDPA_WDOECM_marine.zip",
            gpkg=f"{WDPA_MARINE_DATASET['folder']}/WDPA_WDOECM_marine.gpkg",
        run:
            output_folder = Path(output["zip_file"]).parent
            copy2(input["zip_file"], output["zip_file"])
            unpack_archive(output["zip_file"], output_folder)

            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = f"/vsizip/{output_folder}/WDPA_WDOECM_{bYYYY}_Public_marine_shp_{i}.zip"
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")



if config["enable"]["retrieve"]:

    rule retrieve_monthly_co2_prices:
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
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    rule retrieve_monthly_fuel_prices:
        output:
            "data/validation/energy-price-trends-xlsx-5619002.xlsx",
        log:
            "logs/retrieve_monthly_fuel_prices.log",
        resources:
            mem_mb=5000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_monthly_fuel_prices.py"


if (TYDNP_DATASET := dataset_version("tyndp"))["source"] in ["primary", "archive"]:

    rule retrieve_tyndp:
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
                move(input[key], output[f"{key}_zip"])

                # unzip
                output_folder = Path(output[f"{key}_zip"]).parent
                unpack_archive(output[f"{key}_zip"], output_folder)

                # Remove __MACOSX directory if it exists
                macosx_dir = output_folder / "__MACOSX"
                rmtree(macosx_dir, ignore_errors=True)



if (OSM_DATASET := dataset_version("osm"))["source"] in ["archive"]:
    OSM_FILES = [
        "buses.csv",
        "converters.csv",
        "lines.csv",
        "links.csv",
        "transformers.csv",
        # Newer versions include the additional map.html file for visualisation
        *(["map.html"] if float(OSM_DATASET["version"]) >= 0.6 else []),
    ]

    rule retrieve_osm_archive:
        input:
            **{
                file: storage(f"{OSM_DATASET['url']}/files/{file}")
                for file in OSM_FILES
            },
        output:
            **{file: f"{OSM_DATASET['folder']}/{file}" for file in OSM_FILES},
        log:
            "logs/retrieve_osm_archive.log",
        threads: 1
        resources:
            mem_mb=500,
        run:
            for key in input.keys():
                move(input[key], output[key])
                validate_checksum(output[key], input[key])


elif OSM_DATASET["source"] == "build":

    OSM_FILES = [
        "cables_way.json",
        "lines_way.json",
        "routes_relation.json",
        "substations_way.json",
        "substations_relation.json",
    ]

    rule retrieve_osm_raw:
        output:
            **{
                file.replace(
                    ".json", ""
                ): f"{OSM_DATASET['folder']}/raw/{{country}}/{file}"
                for file in files
            },
        log:
            "logs/retrieve_osm_data_{country}.log",
        threads: 1
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_osm_raw.py"

    rule retrieve_osm_raw_all:
        input:
            expand(
                f"{OSM_DATASET['folder']}/raw/{{country}}/{{file}}",
                country=config_provider("countries"),
                file=OSM_FILES,
            ),


if (NATURA_DATASET := dataset_version("natura"))["source"] in ["archive"]:

    rule retrieve_natura:
        input:
            storage(NATURA_DATASET["url"]),
        output:
            NATURA_DATASET["folder"] / "natura.tiff",
        log:
            "logs/retrieve_natura.log",
        run:
            move(input[0], output[0])

elif NATURA_DATASET["source"] == "build":

    rule build_natura_raster:
        input:
            online=storage(NATURA_DATASET["url"]),
            cutout=lambda w: input_cutout(w),
        output:
            zip=NATURA_DATASET["folder"] / "raw/natura.zip",
            raw=directory(NATURA_DATASET["folder"] / "raw"),
            raster=NATURA_DATASET["folder"] / "natura.tiff",
        resources:
            mem_mb=5000,
        log:
            "logs/build_natura.log",
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_natura.py"



if (OSM_BOUNDARIES_DATASET := dataset_version("osm_boundaries"))["source"] in ["primary"]:

    rule retrieve_osm_boundaries:
        params:
            data_folder=OSM_BOUNDARIES_DATASET["folder"],
            version=OSM_BOUNDARIES_DATASET["version"],
        output:
            json = f"{OSM_BOUNDARIES_DATASET["folder"]}" + "/{country}_adm1.json",
        log:
            "logs/retrieve_osm_boundaries_{country}_adm1.log",
        threads: 1
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_osm_boundaries.py"

elif (OSM_BOUNDARIES_DATASET := dataset_version("osm_boundaries"))["source"] in ["archive"]:
    rule retrieve_osm_boundaries:
        params:
            data_folder=OSM_BOUNDARIES_DATASET["folder"],
            version=OSM_BOUNDARIES_DATASET["version"],
        input:
            storage(
                f"{OSM_BOUNDARIES_DATASET["url"]}",
            ),
        output:
            json1 = f"{OSM_BOUNDARIES_DATASET["folder"]}/XK_adm1.json",
            json2 = f"{OSM_BOUNDARIES_DATASET["folder"]}/UA_adm1.json",
            json3 = f"{OSM_BOUNDARIES_DATASET["folder"]}/MD_adm1.json",
            json4 = f"{OSM_BOUNDARIES_DATASET["folder"]}/BA_adm1.json",
            zip = f"{OSM_BOUNDARIES_DATASET["folder"]}" "/osm_boundaries.zip",
        threads: 1
        run:
            move(input[0],output.zip)
            unpack_archive(output.zip, params.data_folder)


rule retrieve_geothermal_heat_utilisation_potentials:
    input:
        isi_heat_potentials=storage(
            "https://fordatis.fraunhofer.de/bitstream/fordatis/341.5/11/Results_DH_Matching_Cluster.xlsx",
        ),
    output:
        "data/isi_heat_utilisation_potentials.xlsx",
    log:
        "logs/retrieve_geothermal_heat_utilisation_potentials.log",
    threads: 1
    retries: 2
    run:
        move(input[0], output[0])

if (LAU_REGIONS_DATASET := dataset_version("lau_regions"))["source"] in [
    "primary",
    "archive"
]:
    rule retrieve_lau_regions:
        input:
            lau_regions=storage(
                LAU_REGIONS_DATASET["url"]
            ),
        output:
            zip=f"{LAU_REGIONS_DATASET['folder']}/lau_regions.zip",
        log:
            "logs/retrieve_lau_regions.log",
        threads: 1
        retries: 2
        run:
            move(input[0], output[0])



rule retrieve_jrc_ardeco:
    output:
        ardeco_gdp="data/jrc-ardeco/ARDECO-SUVGDP.2021.table.csv",
        ardeco_pop="data/jrc-ardeco/ARDECO-SNPTD.2021.table.csv",
    run:
        urls = {
            "ardeco_gdp": "https://urban.jrc.ec.europa.eu/ardeco-api-v2/rest/export/SUVGDP?version=2021&format=csv-table",
            "ardeco_pop": "https://urban.jrc.ec.europa.eu/ardeco-api-v2/rest/export/SNPTD?version=2021&format=csv-table",
        }

        for key, url in urls.items():
            response = requests.get(url)
            output_path = output[key] if key in urls else None
            if output_path:
                with open(output_path, "wb") as f:
                    f.write(response.content)



if (AQUIFER_DATA_DATASET := dataset_version("aquifer_data"))["source"] in [
        "primary",
        "archive",
    ]:

    rule retrieve_aquifer_data_bgr:
        input:
            zip_file=storage(
                AQUIFER_DATA_DATASET["url"]
            ),
        output:
            aquifer_shapes = expand(f"{AQUIFER_DATA_DATASET['folder']}"+"/ihme1500_aquif_ec4060_v12_poly.{ext}",ext=[
                "shp","shx","dbf","cpg","prj","sbn","sbx",
            ]),
        params:
            files=expand("IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.{ext}", ext=[
                "shp","shx","dbf","cpg","prj","sbn","sbx",
            ])
        run:
            with ZipFile(input.zip_file, "r") as zip_ref:
                for fn, outpt in zip(
                    params["files"],
                    output["aquifer_shapes"],
                ):  
                    zip_ref.extract(fn, Path(outpt).parent)
                    extracted_file = Path(outpt).parent / fn
                    extracted_file.rename(outpt)


if config["enable"]["retrieve"]:

    rule retrieve_dh_areas:
        input:
            dh_areas=storage(
                "https://fordatis.fraunhofer.de/bitstream/fordatis/341.5/2/dh_areas.gpkg",
            ),
        output:
            dh_areas="data/dh_areas.gpkg",
        log:
            "logs/retrieve_dh_areas.log",
        threads: 1
        retries: 2
        run:
            move(input[0], output[0])


if (MOBILITY_PROFILES_DATASET := dataset_version("mobility_profiles"))["source"] in [
    "archive"
]:

    rule retrieve_mobility_profiles:
        input:
            kfz=storage(MOBILITY_PROFILES_DATASET["url"] + "/kfz.csv"),
            pkw=storage(MOBILITY_PROFILES_DATASET["url"] + "/pkw.csv"),
        output:
            kfz=MOBILITY_PROFILES_DATASET["folder"] / "kfz.csv",
            pkw=MOBILITY_PROFILES_DATASET["folder"] / "pkw.csv",
        threads: 1
        resources:
            mem_mb=1000,
        log:
            "logs/retrieve_mobility_profiles.log",
        benchmark:
            "benchmarks/retrieve_mobility_profiles"
        conda:
            "../envs/environment.yaml"
        run:
            move(input["kfz"], output["kfz"])
            move(input["pkw"], output["pkw"])
