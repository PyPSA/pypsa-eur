# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import requests
from datetime import datetime, timedelta

if config["enable"].get("retrieve", "auto") == "auto":
    config["enable"]["retrieve"] = has_internet_access()

if config["enable"]["retrieve"] is False:
    print("Datafile downloads disabled in config[retrieve] or no internet access.")


if config["enable"]["retrieve"] and config["enable"].get("retrieve_databundle", True):
    datafiles = [
        "je-e-21.03.02.xls",
        "eez/World_EEZ_v8_2014.shp",
        "naturalearth/ne_10m_admin_0_countries.shp",
        "NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp",
        "nama_10r_3popgdp.tsv.gz",
        "nama_10r_3gdp.tsv.gz",
        "corine/g250_clc06_V18_5.tif",
        "eea/UNFCCC_v23.csv",
        "nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson",
        "myb1-2017-nitro.xls",
        "emobility/KFZ__count",
        "emobility/Pkw__count",
        "h2_salt_caverns_GWh_per_sqkm.geojson",
        "natura/natura.tiff",
        "gebco/GEBCO_2014_2D.nc",
    ]

    rule retrieve_databundle:
        output:
            expand("data/bundle/{file}", file=datafiles),
            directory("data/bundle/jrc-idees-2015"),
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
        script:
            "../scripts/retrieve_eurostat_data.py"

    rule retrieve_eurostat_household_data:
        output:
            "data/eurostat/eurostat-household_energy_balances-february_2024.csv",
        log:
            "logs/retrieve_eurostat_household_data.log",
        retries: 2
        script:
            "../scripts/retrieve_eurostat_household_data.py"


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            storage(
                "https://zenodo.org/records/6382570/files/{cutout}.nc",
            ),
        output:
            protected("cutouts/" + CDIR + "{cutout}.nc"),
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
                "https://zenodo.org/records/10973944/files/shipdensity_global.zip",
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
