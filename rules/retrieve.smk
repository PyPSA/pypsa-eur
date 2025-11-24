# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import os
import requests
from datetime import datetime, timedelta
from shutil import move, unpack_archive
from shutil import copy2 as shcopy2
from zipfile import ZipFile


# Configure the default storage provider for accessing remote files using http
# and the special storage plugin for accessing Zenodo files
storage:
    provider="http",
    retries=3,


storage cached_http:
    provider="cached-http",


if config["enable"].get("retrieve", "auto") == "auto":
    config["enable"]["retrieve"] = has_internet_access()

if config["enable"]["retrieve"] is False:
    print("Datafile downloads disabled in config[retrieve] or no internet access.")


if config["enable"]["retrieve"] and config["enable"].get("retrieve_databundle", True):
    datafiles = [
        "je-e-21.03.02.xls",
        "nama_10r_3popgdp.tsv.gz",
        "corine/g250_clc06_V18_5.tif",
        "eea/UNFCCC_v23.csv",
        "emobility/KFZ__count",
        "emobility/Pkw__count",
        "h2_salt_caverns_GWh_per_sqkm.geojson",
        "natura/natura.tiff",
        "gebco/GEBCO_2014_2D.nc",
        "GDP_per_capita_PPP_1990_2015_v2.nc",
        "ppp_2019_1km_Aggregated.tif",
        "era5-HDD-per-country.csv",
        "era5-runoff-per-country.csv",
    ]

    rule retrieve_databundle:
        output:
            expand("data/bundle/{file}", file=datafiles),
        log:
            "logs/retrieve_databundle.log",
        benchmark:
            "benchmarks/retrieve_databundle"
        resources:
            mem_mb=1000,
        retries: 2
        script:
            "../scripts/retrieve_databundle.py"

    rule retrieve_eurostat_data:
        output:
            directory("data/eurostat/Balances-April2023"),
        log:
            "logs/retrieve_eurostat_data.log",
        retries: 2
        script:
            "../scripts/retrieve_eurostat_data.py"

    rule retrieve_jrc_idees:
        output:
            directory("data/jrc-idees-2021"),
        log:
            "logs/retrieve_jrc_idees.log",
        retries: 2
        script:
            "../scripts/retrieve_jrc_idees.py"

    rule retrieve_eurostat_household_data:
        output:
            "data/eurostat/eurostat-household_energy_balances-february_2024.csv",
        log:
            "logs/retrieve_eurostat_household_data.log",
        retries: 2
        script:
            "../scripts/retrieve_eurostat_household_data.py"


if config["enable"]["retrieve"]:

    rule retrieve_nuts_2013_shapes:
        input:
            shapes=storage(
                "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2013-03m.geojson.zip"
            ),
        output:
            shapes_level_3="data/nuts/NUTS_RG_03M_2013_4326_LEVL_3.geojson",
            shapes_level_2="data/nuts/NUTS_RG_03M_2013_4326_LEVL_2.geojson",
        params:
            zip_file="ref-nuts-2013-03m.geojson.zip",
        run:
            # Copy file and ensure proper permissions
            shcopy2(input.shapes, params.zip_file)

            with ZipFile(params.zip_file, "r") as zip_ref:
                for level in ["LEVL_3", "LEVL_2"]:
                    filename = f"NUTS_RG_03M_2013_4326_{level}.geojson"
                    zip_ref.extract(filename, Path(output.shapes_level_3).parent)
                    extracted_file = Path(output.shapes_level_3).parent / filename
                    extracted_file.rename(
                        getattr(output, f"shapes_level_{level[-1]}")
                    )
            os.remove(params.zip_file)



if config["enable"]["retrieve"]:

    rule retrieve_nuts_2021_shapes:
        input:
            shapes=storage(
                "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2021-01m.geojson.zip"
            ),
        output:
            shapes_level_3="data/nuts/NUTS_RG_01M_2021_4326_LEVL_3.geojson",
            shapes_level_2="data/nuts/NUTS_RG_01M_2021_4326_LEVL_2.geojson",
            shapes_level_1="data/nuts/NUTS_RG_01M_2021_4326_LEVL_1.geojson",
            shapes_level_0="data/nuts/NUTS_RG_01M_2021_4326_LEVL_0.geojson",
        params:
            zip_file="ref-nuts-2021-01m.geojson.zip",
        run:
            # Copy file and ensure proper permissions
            shcopy2(input.shapes, params.zip_file)

            with ZipFile(params.zip_file, "r") as zip_ref:
                for level in ["LEVL_3", "LEVL_2", "LEVL_1", "LEVL_0"]:
                    filename = f"NUTS_RG_01M_2021_4326_{level}.geojson"
                    zip_ref.extract(filename, Path(output.shapes_level_0).parent)
                    extracted_file = Path(output.shapes_level_0).parent / filename
                    extracted_file.rename(
                        getattr(output, f"shapes_level_{level[-1]}")
                    )
            os.remove(params.zip_file)



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
        script:
            "../scripts/retrieve_bidding_zones.py"


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            storage(
                "https://zenodo.org/records/15349674/files/{cutout}.nc",
            ),
        output:
            CDIR.joinpath("{cutout}.nc").as_posix(),
        log:
            Path("logs").joinpath(CDIR, "retrieve_cutout_{cutout}.log").as_posix(),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])
            validate_checksum(output[0], input[0])


if config["enable"]["retrieve"]:

    rule retrieve_tyndp_bundle:
        output:
            reference_grid="data/tyndp_2024_bundle/Line data/ReferenceGrid_Electricity.xlsx",
            buses="data/tyndp_2024_bundle/Nodes/LIST OF NODES.xlsx",
        log:
            "logs/retrieve_tyndp_bundle.log",
        retries: 2
        script:
            "../scripts/retrieve_tyndp_bundle.py"


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        params:
            version=config_provider("costs", "version"),
        output:
            resources("costs_{planning_horizons}.csv"),
        log:
            logs("retrieve_cost_data_{planning_horizons}.log"),
        resources:
            mem_mb=1000,
        retries: 2
        script:
            "../scripts/retrieve_cost_data.py"


if config["enable"]["retrieve"]:
    datafiles = [
        "IGGIELGN_LNGs.geojson",
        "IGGIELGN_BorderPoints.geojson",
        "IGGIELGN_Productions.geojson",
        "IGGIELGN_Storages.geojson",
        "IGGIELGN_PipeSegments.geojson",
    ]

    rule retrieve_gas_infrastructure_data:
        output:
            expand("data/gas_network/scigrid-gas/data/{files}", files=datafiles),
        log:
            "logs/retrieve_gas_infrastructure_data.log",
        retries: 2
        script:
            "../scripts/retrieve_gas_infrastructure_data.py"


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
        script:
            "../scripts/retrieve_electricity_demand.py"


if config["enable"]["retrieve"]:

    rule retrieve_synthetic_electricity_demand:
        input:
            storage(
                "https://zenodo.org/records/10820928/files/demand_hourly.csv",
            ),
        output:
            "data/load_synthetic_raw.csv",
        log:
            "logs/retrieve_synthetic_electricity_demand.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    rule retrieve_ship_raster:
        input:
            storage(
                "https://zenodo.org/records/13757228/files/shipdensity_global.zip",
                keep_local=True,
            ),
        output:
            "data/shipdensity_global.zip",
        log:
            "logs/retrieve_ship_raster.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])
            validate_checksum(output[0], input[0])


if config["enable"]["retrieve"]:

    rule retrieve_jrc_enspreso_biomass:
        input:
            storage(
                "https://zenodo.org/records/10356004/files/ENSPRESO_BIOMASS.xlsx",
                keep_local=True,
            ),
        output:
            "data/ENSPRESO_BIOMASS.xlsx",
        retries: 1
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    rule retrieve_hotmaps_industrial_sites:
        input:
            storage(
                "https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database/-/raw/master/data/Industrial_Database.csv",
                keep_local=True,
            ),
        output:
            "data/Industrial_Database.csv",
        retries: 1
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    rule retrieve_usgs_ammonia_production:
        input:
            storage(
                "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/myb1-2022-nitro-ert.xlsx"
            ),
        output:
            "data/myb1-2022-nitro-ert.xlsx",
        retries: 1
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    # Downloading Copernicus Global Land Cover for land cover and land use:
    # Website: https://land.copernicus.eu/global/products/lc
    rule download_copernicus_land_cover:
        input:
            storage(
                "https://zenodo.org/records/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
            ),
        output:
            "data/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        run:
            move(input[0], output[0])
            validate_checksum(output[0], input[0])


if config["enable"]["retrieve"]:

    # Downloading LUISA Base Map for land cover and land use:
    # Website: https://ec.europa.eu/jrc/en/luisa
    rule retrieve_luisa_land_cover:
        input:
            storage(
                "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/LUISA/EUROPE/Basemaps/LandUse/2018/LATEST/LUISA_basemap_020321_50m.tif",
            ),
        output:
            "data/LUISA_basemap_020321_50m.tif",
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    rule retrieve_eez:
        params:
            zip_file="World_EEZ_v12_20231025_LR.zip",
        output:
            gpkg="data/eez/World_EEZ_v12_20231025_LR/eez_v12_lowres.gpkg",
        run:
            from uuid import uuid4

            name = str(uuid4())[:8]
            org = str(uuid4())[:8]

            response = requests.post(
                "https://www.marineregions.org/download_file.php",
                params={"name": "World_EEZ_v12_20231025_LR.zip"},
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

            with open(params["zip_file"], "wb") as f:
                f.write(response.content)
            output_folder = Path(output.gpkg).parent.parent
            unpack_archive(params["zip_file"], output_folder)
            os.remove(params["zip_file"])



if config["enable"]["retrieve"]:

    rule retrieve_worldbank_urban_population:
        params:
            zip_file="API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.zip",
        output:
            gpkg="data/worldbank/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.csv",
        run:
            response = requests.get(
                "https://api.worldbank.org/v2/en/indicator/SP.URB.TOTL.IN.ZS?downloadformat=csv",
            )

            with open(params["zip_file"], "wb") as f:
                f.write(response.content)
            output_folder = Path(output.gpkg).parent
            unpack_archive(params["zip_file"], output_folder)

            for f in os.listdir(output_folder):
                if f.startswith(
                    "API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2_"
                ) and f.endswith(".csv"):
                    os.rename(os.path.join(output_folder, f), output.gpkg)
                    break
            os.remove(params["zip_file"])



if config["enable"]["retrieve"]:

    rule retrieve_co2stop:
        params:
            zip_file="co2jrc_openformats.zip",
        output:
            "data/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Storage_Units.csv",
            "data/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps.csv",
            "data/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps_Temp.csv",
            "data/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps1.csv",
            "data/CO2JRC_OpenFormats/CO2Stop_Polygons Data/DaughterUnits_March13.kml",
            "data/CO2JRC_OpenFormats/CO2Stop_Polygons Data/StorageUnits_March13.kml",
        run:
            response = requests.get(
                "https://setis.ec.europa.eu/document/download/786a884f-0b33-4789-b744-28004b16bd1a_en?filename=co2jrc_openformats.zip",
            )
            with open(params["zip_file"], "wb") as f:
                f.write(response.content)
            output_folder = Path(output[0]).parent.parent.parent
            unpack_archive(params["zip_file"], output_folder)
            os.remove(params["zip_file"])



if config["enable"]["retrieve"]:

    rule retrieve_gem_europe_gas_tracker:
        output:
            "data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        run:
            # mirror of https://globalenergymonitor.org/wp-content/uploads/2024/05/Europe-Gas-Tracker-2024-05.xlsx
            url = "https://tubcloud.tu-berlin.de/s/LMBJQCsN6Ez5cN2/download/Europe-Gas-Tracker-2024-05.xlsx"
            response = requests.get(url)
            with open(output[0], "wb") as f:
                f.write(response.content)



if config["enable"]["retrieve"]:

    rule retrieve_gem_steel_plant_tracker:
        output:
            "data/gem/Global-Steel-Plant-Tracker-April-2024-Standard-Copy-V1.xlsx",
        run:
            # mirror or https://globalenergymonitor.org/wp-content/uploads/2024/04/Global-Steel-Plant-Tracker-April-2024-Standard-Copy-V1.xlsx
            url = "https://tubcloud.tu-berlin.de/s/Aqebo3rrQZWKGsG/download/Global-Steel-Plant-Tracker-April-2024-Standard-Copy-V1.xlsx"
            response = requests.get(url)
            with open(output[0], "wb") as f:
                f.write(response.content)



if config["enable"]["retrieve"]:
    # Some logic to find the correct file URL
    # Sometimes files are released delayed or ahead of schedule, check which file is currently available

    def check_file_exists(url):
        response = requests.head(url)
        return response.status_code == 200

    # Basic pattern where WDPA files can be found
    url_pattern = (
        "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_{bYYYY}_Public_shp.zip"
    )

    # 3-letter month + 4 digit year for current/previous/next month to test
    current_monthyear = datetime.now().strftime("%b%Y")
    prev_monthyear = (datetime.now() - timedelta(30)).strftime("%b%Y")
    next_monthyear = (datetime.now() + timedelta(30)).strftime("%b%Y")

    # Test prioritised: current month -> previous -> next
    for bYYYY in [current_monthyear, prev_monthyear, next_monthyear]:
        if check_file_exists(url := url_pattern.format(bYYYY=bYYYY)):
            break
        else:
            # If None of the three URLs are working
            url = False

    # assert (
    #     url
    # ), f"No WDPA files found at {url_pattern} for bY='{current_monthyear}, {prev_monthyear}, or {next_monthyear}'"

    # Downloading protected area database from WDPA
    # extract the main zip and then merge the contained 3 zipped shapefiles
    # Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
    rule download_wdpa:
        input:
            zip_file=storage(url, keep_local=True) if url else [],
        params:
            zip_file="WDPA_shp.zip",
            folder_name="WDPA",
        output:
            gpkg="data/WDPA.gpkg",
        run:
            # Copy file and ensure proper permissions
            shcopy2(input.zip_file, params.zip_file)
            output_folder = Path(output.gpkg).parent / params.folder_name
            unpack_archive(params.zip_file, output_folder)

            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = (
                    f"/vsizip/{output_folder}/WDPA_{bYYYY}_Public_shp_{i}.zip"
                )
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")
            os.remove(params.zip_file)

    rule download_wdpa_marine:
        # Downloading Marine protected area database from WDPA
        # extract the main zip and then merge the contained 3 zipped shapefiles
        # Website: https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
        input:
            zip_file=storage(
                f"https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{bYYYY}_Public_marine_shp.zip",
                keep_local=True,
            ),
        params:
            zip_file="WDPA_WDOECM_marine.zip",
            folder_name="WDPA_WDOECM_marine",
        output:
            gpkg="data/WDPA_WDOECM_marine.gpkg",
        run:
            shcopy2(input.zip_file, params.zip_file)
            output_folder = Path(output.gpkg).parent / params.folder_name
            unpack_archive(params.zip_file, output_folder)

            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = f"/vsizip/{output_folder}/WDPA_WDOECM_{bYYYY}_Public_marine_shp_{i}.zip"
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")
            os.remove(params.zip_file)



if config["enable"]["retrieve"]:

    rule retrieve_monthly_co2_prices:
        input:
            storage(
                "https://public.eex-group.com/eex/eua-auction-report/emission-spot-primary-market-auction-report-2019-data.xls",
                keep_local=True,
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
        script:
            "../scripts/retrieve_monthly_fuel_prices.py"


if config["enable"]["retrieve"] and (
    config["electricity"]["base_network"] == "osm-prebuilt"
):
    OSM_VERSION = config["electricity"]["osm-prebuilt-version"]
    OSM_FILES = [
        "buses.csv",
        "converters.csv",
        "lines.csv",
        "links.csv",
        "transformers.csv",
    ]
    if OSM_VERSION >= 0.6:
        OSM_FILES.append("map.html")
    OSM_ZENODO_IDS = {
        0.1: "12799202",
        0.2: "13342577",
        0.3: "13358976",
        0.4: "13759222",
        0.5: "13981528",
        0.6: "14144752",
    }

    # update rule to use the correct version
    rule retrieve_osm_prebuilt:
        input:
            **{
                file: storage(
                    f"https://zenodo.org/records/{OSM_ZENODO_IDS[OSM_VERSION]}/files/{file}"
                )
                for file in OSM_FILES
            },
        output:
            **{file: f"data/osm-prebuilt/{OSM_VERSION}/{file}" for file in OSM_FILES},
        log:
            "logs/retrieve_osm_prebuilt.log",
        threads: 1
        resources:
            mem_mb=500,
        run:
            for key in input.keys():
                move(input[key], output[key])
                validate_checksum(output[key], input[key])



if config["enable"]["retrieve"] and (
    config["electricity"]["base_network"] == "osm-raw"
):

    rule retrieve_osm_data:
        output:
            cables_way="data/osm-raw/{country}/cables_way.json",
            lines_way="data/osm-raw/{country}/lines_way.json",
            routes_relation="data/osm-raw/{country}/routes_relation.json",
            substations_way="data/osm-raw/{country}/substations_way.json",
            substations_relation="data/osm-raw/{country}/substations_relation.json",
        log:
            "logs/retrieve_osm_data_{country}.log",
        threads: 1
        script:
            "../scripts/retrieve_osm_data.py"


if config["enable"]["retrieve"] and (
    config["electricity"]["base_network"] == "osm-raw"
):

    rule retrieve_osm_data_all:
        input:
            expand(
                "data/osm-raw/{country}/cables_way.json",
                country=config_provider("countries"),
            ),
            expand(
                "data/osm-raw/{country}/lines_way.json",
                country=config_provider("countries"),
            ),
            expand(
                "data/osm-raw/{country}/routes_relation.json",
                country=config_provider("countries"),
            ),
            expand(
                "data/osm-raw/{country}/substations_way.json",
                country=config_provider("countries"),
            ),
            expand(
                "data/osm-raw/{country}/substations_relation.json",
                country=config_provider("countries"),
            ),


if config["enable"]["retrieve"]:

    rule retrieve_osm_boundaries:
        output:
            json="data/osm-boundaries/json/{country}_adm1.json",
        log:
            "logs/retrieve_osm_boundaries_{country}_adm1.log",
        threads: 1
        script:
            "../scripts/retrieve_osm_boundaries.py"

    rule retrieve_geothermal_heat_utilisation_potentials:
        input:
            isi_heat_potentials=storage(
                "https://fordatis.fraunhofer.de/bitstream/fordatis/341.5/11/Results_DH_Matching_Cluster.xlsx",
                keep_local=True,
            ),
        output:
            "data/isi_heat_utilisation_potentials.xlsx",
        log:
            "logs/retrieve_geothermal_heat_utilisation_potentials.log",
        threads: 1
        retries: 2
        run:
            move(input[0], output[0])

    rule retrieve_lau_regions:
        input:
            lau_regions=storage(
                "https://gisco-services.ec.europa.eu/distribution/v2/lau/download/ref-lau-2019-01m.geojson.zip",
                keep_local=True,
            ),
        output:
            lau_regions="data/lau_regions.zip",
        log:
            "logs/retrieve_lau_regions.log",
        log:
            "logs/retrieve_lau_regions.log",
        threads: 1
        retries: 2
        run:
            move(input[0], output[0])

    rule retrieve_seawater_temperature:
        params:
            default_cutout=config_provider("atlite", "default_cutout"),
        output:
            seawater_temperature="data/seawater_temperature_{year}.nc",
        log:
            "logs/retrieve_seawater_temperature_{year}.log",
        resources:
            mem_mb=10000,
        script:
            "../scripts/retrieve_seawater_temperature.py"

    rule retrieve_hera_data_test_cutout:
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


if config["enable"]["retrieve"]:

    rule retrieve_jrc_ardeco:
        output:
            ardeco_gdp="data/jrc-ardeco/ARDECO-SUVGDP.2021.table.csv",
            ardeco_pop="data/jrc-ardeco/ARDECO-SNPTD.2021.table.csv",
        run:
            urls = {
                "ardeco_gdp": "https://territorial.ec.europa.eu/ardeco-api-v2/rest/export/SUVGDP?versions=2021&unit=EUR&level_id=0&level_id=1&level_id=2&level_id=3&format=csv-table",
                "ardeco_pop": "https://territorial.ec.europa.eu/ardeco-api-v2/rest/export/SNPTD?versions=2021&unit=NR&level_id=0&level_id=1&level_id=2&level_id=3&format=csv-table",
            }

            for key, url in urls.items():
                response = requests.get(url)
                output_path = output[key] if key in urls else None
                if output_path:
                    with open(output_path, "wb") as f:
                        f.write(response.content)



if config["enable"]["retrieve"]:

    rule retrieve_aquifer_data_bgr:
        input:
            zip_file=storage(
                "https://download.bgr.de/bgr/grundwasser/IHME1500/v12/shp/IHME1500_v12.zip"
            ),
        output:
            aquifer_shapes_shp="data/bgr/ihme1500_aquif_ec4060_v12_poly.shp",
            aquifer_shapes_shx="data/bgr/ihme1500_aquif_ec4060_v12_poly.shx",
            aquifer_shapes_dbf="data/bgr/ihme1500_aquif_ec4060_v12_poly.dbf",
            aquifer_shapes_cpg="data/bgr/ihme1500_aquif_ec4060_v12_poly.cpg",
            aquifer_shapes_prj="data/bgr/ihme1500_aquif_ec4060_v12_poly.prj",
            aquifer_shapes_sbn="data/bgr/ihme1500_aquif_ec4060_v12_poly.sbn",
            aquifer_shapes_sbx="data/bgr/ihme1500_aquif_ec4060_v12_poly.sbx",
        params:
            filename_shp="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.shp",
            filename_shx="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.shx",
            filename_dbf="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.dbf",
            filename_cpg="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.cpg",
            filename_prj="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.prj",
            filename_sbn="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.sbn",
            filename_sbx="IHME1500_v12/shp/ihme1500_aquif_ec4060_v12_poly.sbx",
        run:
            with ZipFile(input.zip_file, "r") as zip_ref:
                for fn, outpt in zip(
                    params,
                    output,
                ):
                    zip_ref.extract(fn, Path(outpt).parent)
                    extracted_file = Path(outpt).parent / fn
                    extracted_file.rename(outpt)

    rule retrieve_dh_areas:
        input:
            dh_areas=storage(
                "https://fordatis.fraunhofer.de/bitstream/fordatis/341.5/2/dh_areas.gpkg",
                keep_local=True,
            ),
        output:
            dh_areas="data/dh_areas.gpkg",
        log:
            "logs/retrieve_dh_areas.log",
        threads: 1
        retries: 2
        run:
            move(input[0], output[0])
