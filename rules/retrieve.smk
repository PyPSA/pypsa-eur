# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import requests
from datetime import datetime, timedelta
from shutil import move, unpack_archive
from zipfile import ZipFile

if config["enable"].get("retrieve", "auto") == "auto":
    config["enable"]["retrieve"] = has_internet_access()

if config["enable"]["retrieve"] is False:
    print("Datafile downloads disabled in config[retrieve] or no internet access.")


if config["enable"]["retrieve"] and config["enable"].get("retrieve_databundle", True):
    datafiles = [
        "je-e-21.03.02.xls",
        "nama_10r_3popgdp.tsv.gz",
        "nama_10r_3gdp.tsv.gz",
        "corine/g250_clc06_V18_5.tif",
        "eea/UNFCCC_v23.csv",
        "emobility/KFZ__count",
        "emobility/Pkw__count",
        "h2_salt_caverns_GWh_per_sqkm.geojson",
        "natura/natura.tiff",
        "gebco/GEBCO_2014_2D.nc",
        "GDP_per_capita_PPP_1990_2015_v2.nc",
        "ppp_2013_1km_Aggregated.tif",
    ]

    rule retrieve_databundle:
        output:
            expand("data/bundle/{file}", file=datafiles),
        log:
            "logs/retrieve_databundle.log",
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/retrieve.yaml"
        script:
            "../scripts/retrieve_databundle.py"

    rule retrieve_eurostat_data:
        output:
            directory("data/eurostat/Balances-April2023"),
        log:
            "logs/retrieve_eurostat_data.log",
        retries: 2
        conda:
            "../envs/retrieve.yaml"
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
        conda:
            "../envs/retrieve.yaml"
        script:
            "../scripts/retrieve_eurostat_household_data.py"


if config["enable"]["retrieve"]:

    rule retrieve_nuts_shapes:
        input:
            shapes=storage(
                "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2013-03m.geojson.zip"
            ),
        output:
            shapes_level_3="data/nuts/NUTS_RG_03M_2013_4326_LEVL_3.geojson",
            shapes_level_2="data/nuts/NUTS_RG_03M_2013_4326_LEVL_2.geojson",
        params:
            zip_file="data/nuts/ref-nuts-2013-03m.geojson.zip",
        run:
            os.rename(input.shapes, params.zip_file)
            with ZipFile(params.zip_file, "r") as zip_ref:
                for level in ["LEVL_3", "LEVL_2"]:
                    filename = f"NUTS_RG_03M_2013_4326_{level}.geojson"
                    zip_ref.extract(filename, Path(output.shapes_level_3).parent)
                    extracted_file = Path(output.shapes_level_3).parent / filename
                    extracted_file.rename(
                        getattr(output, f"shapes_level_{level[-1]}")
                    )
            os.remove(params.zip_file)



if config["enable"]["retrieve"] and config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            storage(
                "https://zenodo.org/records/12791128/files/{cutout}.nc",
            ),
        output:
            "cutouts/" + CDIR + "{cutout}.nc",
        log:
            "logs/" + CDIR + "retrieve_cutout_{cutout}.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])
            validate_checksum(output[0], input[0])


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        params:
            version=config_provider("costs", "version"),
        output:
            resources("costs_{year}.csv"),
        log:
            logs("retrieve_cost_data_{year}.log"),
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/retrieve.yaml"
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
        conda:
            "../envs/retrieve.yaml"
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
        conda:
            "../envs/retrieve.yaml"
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

    rule retrieve_geological_co2_storage_potential:
        input:
            storage(
                "https://raw.githubusercontent.com/ericzhou571/Co2Storage/main/resources/complete_map_2020_unit_Mt.geojson",
                keep_local=True,
            ),
        output:
            "data/complete_map_2020_unit_Mt.geojson",
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
            zip="data/eez/World_EEZ_v12_20231025_LR.zip",
        output:
            gpkg="data/eez/World_EEZ_v12_20231025_LR/eez_v12_lowres.gpkg",
        run:
            import os
            import requests
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

            with open(params["zip"], "wb") as f:
                f.write(response.content)
            output_folder = Path(params["zip"]).parent
            unpack_archive(params["zip"], output_folder)
            os.remove(params["zip"])



if config["enable"]["retrieve"]:

    rule retrieve_worldbank_urban_population:
        params:
            zip="data/worldbank/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.zip",
        output:
            gpkg="data/worldbank/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.csv",
        run:
            import os
            import requests

            response = requests.get(
                "https://api.worldbank.org/v2/en/indicator/SP.URB.TOTL.IN.ZS?downloadformat=csv",
            )

            with open(params["zip"], "wb") as f:
                f.write(response.content)
            output_folder = Path(params["zip"]).parent
            unpack_archive(params["zip"], output_folder)

            for f in os.listdir(output_folder):
                if f.startswith(
                    "API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2_"
                ) and f.endswith(".csv"):
                    os.rename(os.path.join(output_folder, f), output.gpkg)
                    break



if config["enable"]["retrieve"]:

    # Download directly from naciscdn.org which is a redirect from naturalearth.com
    # (https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/)
    # Use point-of-view (POV) variant of Germany so that Crimea is included.
    rule retrieve_naturalearth_countries:
        input:
            storage(
                "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries_deu.zip"
            ),
        params:
            zip="data/naturalearth/ne_10m_admin_0_countries_deu.zip",
        output:
            countries="data/naturalearth/ne_10m_admin_0_countries_deu.shp",
        run:
            move(input[0], params["zip"])
            output_folder = Path(output["countries"]).parent
            unpack_archive(params["zip"], output_folder)
            os.remove(params["zip"])


if config["enable"]["retrieve"]:

    rule retrieve_gem_europe_gas_tracker:
        output:
            "data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        run:
            import requests

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
            import requests

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

    assert (
        url
    ), f"No WDPA files found at {url_pattern} for bY='{current_monthyear}, {prev_monthyear}, or {next_monthyear}'"

    # Downloading protected area database from WDPA
    # extract the main zip and then merge the contained 3 zipped shapefiles
    # Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
    rule download_wdpa:
        input:
            storage(url, keep_local=True),
        params:
            zip="data/WDPA_shp.zip",
            folder=directory("data/WDPA"),
        output:
            gpkg="data/WDPA.gpkg",
        run:
            shell("cp {input} {params.zip}")
            shell("unzip -o {params.zip} -d {params.folder}")
            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = (
                    f"/vsizip/{params.folder}/WDPA_{bYYYY}_Public_shp_{i}.zip"
                )
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")

    rule download_wdpa_marine:
        # Downloading Marine protected area database from WDPA
        # extract the main zip and then merge the contained 3 zipped shapefiles
        # Website: https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
        input:
            storage(
                f"https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{bYYYY}_Public_marine_shp.zip",
                keep_local=True,
            ),
        params:
            zip="data/WDPA_WDOECM_marine.zip",
            folder=directory("data/WDPA_WDOECM_marine"),
        output:
            gpkg="data/WDPA_WDOECM_marine.gpkg",
        run:
            shell("cp {input} {params.zip}")
            shell("unzip -o {params.zip} -d {params.folder}")
            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = f"/vsizip/{params.folder}/WDPA_WDOECM_{bYYYY}_Public_marine_shp_{i}.zip"
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")



if config["enable"]["retrieve"]:

    rule retrieve_monthly_co2_prices:
        input:
            storage(
                "https://www.eex.com/fileadmin/EEX/Downloads/EUA_Emission_Spot_Primary_Market_Auction_Report/Archive_Reports/emission-spot-primary-market-auction-report-2019-data.xls",
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
        conda:
            "../envs/retrieve.yaml"
        script:
            "../scripts/retrieve_monthly_fuel_prices.py"


if config["enable"]["retrieve"] and (
    config["electricity"]["base_network"] == "osm-prebuilt"
):
    # Dictionary of prebuilt versions, e.g. 0.3 : "13358976"
    osm_prebuilt_version = {
        0.1: "12799202",
        0.2: "13342577",
        0.3: "13358976",
        0.4: "13759222",
        0.5: "13981528",
    }

    # update rule to use the correct version
    rule retrieve_osm_prebuilt:
        input:
            buses=storage(
                f"https://zenodo.org/records/{osm_prebuilt_version[config['electricity']['osm-prebuilt-version']]}/files/buses.csv"
            ),
            converters=storage(
                f"https://zenodo.org/records/{osm_prebuilt_version[config['electricity']['osm-prebuilt-version']]}/files/converters.csv"
            ),
            lines=storage(
                f"https://zenodo.org/records/{osm_prebuilt_version[config['electricity']['osm-prebuilt-version']]}/files/lines.csv"
            ),
            links=storage(
                f"https://zenodo.org/records/{osm_prebuilt_version[config['electricity']['osm-prebuilt-version']]}/files/links.csv"
            ),
            transformers=storage(
                f"https://zenodo.org/records/{osm_prebuilt_version[config['electricity']['osm-prebuilt-version']]}/files/transformers.csv"
            ),
        output:
            buses=f"data/osm-prebuilt/{config['electricity']['osm-prebuilt-version']}/buses.csv",
            converters=f"data/osm-prebuilt/{config['electricity']['osm-prebuilt-version']}/converters.csv",
            lines=f"data/osm-prebuilt/{config['electricity']['osm-prebuilt-version']}/lines.csv",
            links=f"data/osm-prebuilt/{config['electricity']['osm-prebuilt-version']}/links.csv",
            transformers=f"data/osm-prebuilt/{config['electricity']['osm-prebuilt-version']}/transformers.csv",
        log:
            "logs/retrieve_osm_prebuilt.log",
        threads: 1
        resources:
            mem_mb=500,
        retries: 2
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
        conda:
            "../envs/retrieve.yaml"
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
