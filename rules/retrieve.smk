# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
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
        "ch_cantons.csv",
        "je-e-21.03.02.xls",
        "eez/World_EEZ_v8_2014.shp",
        "hydro_capacities.csv",
        "naturalearth/ne_10m_admin_0_countries.shp",
        "NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp",
        "nama_10r_3popgdp.tsv.gz",
        "nama_10r_3gdp.tsv.gz",
        "corine/g250_clc06_V18_5.tif",
    ]

    if not config.get("tutorial", False):
        datafiles.extend(["natura/Natura2000_end2015.shp", "GEBCO_2014_2D.nc"])

    rule retrieve_databundle:
        output:
            protected(expand("data/bundle/{file}", file=datafiles)),
        log:
            LOGS + "retrieve_databundle.log",
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_databundle.py"


if config["enable"].get("retrieve_irena"):

    rule retrieve_irena:
        output:
            offwind="data/existing_infrastructure/offwind_capacity_IRENA.csv",
            onwind="data/existing_infrastructure/onwind_capacity_IRENA.csv",
            solar="data/existing_infrastructure/solar_capacity_IRENA.csv",
        log:
            LOGS + "retrieve_irena.log",
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_irena.py"


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            HTTP.remote(
                "zenodo.org/record/6382570/files/{cutout}.nc",
                static=True,
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


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        input:
            HTTP.remote(
                "raw.githubusercontent.com/PyPSA/technology-data/{}/outputs/".format(
                    config["costs"]["version"]
                )
                + "costs_{year}.csv",
                keep_local=True,
            ),
        output:
            "data/costs_{year}.csv",
        log:
            LOGS + "retrieve_cost_data_{year}.log",
        resources:
            mem_mb=1000,
        retries: 2
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"] and config["enable"].get(
    "retrieve_natura_raster", True
):

    rule retrieve_natura_raster:
        input:
            HTTP.remote(
                "zenodo.org/record/4706686/files/natura.tiff",
                keep_local=True,
                static=True,
            ),
        output:
            protected(RESOURCES + "natura.tiff"),
        log:
            LOGS + "retrieve_natura_raster.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"] and config["enable"].get(
    "retrieve_sector_databundle", True
):
    datafiles = [
        "eea/UNFCCC_v23.csv",
        "switzerland-sfoe/switzerland-new_format.csv",
        "nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson",
        "myb1-2017-nitro.xls",
        "Industrial_Database.csv",
        "emobility/KFZ__count",
        "emobility/Pkw__count",
        "h2_salt_caverns_GWh_per_sqkm.geojson",
    ]

    datafolders = [
        protected(
            directory("data/bundle-sector/eurostat-energy_balances-june_2016_edition")
        ),
        protected(
            directory("data/bundle-sector/eurostat-energy_balances-may_2018_edition")
        ),
        protected(directory("data/bundle-sector/jrc-idees-2015")),
    ]

    rule retrieve_sector_databundle:
        output:
            protected(expand("data/bundle-sector/{files}", files=datafiles)),
            *datafolders,
        log:
            LOGS + "retrieve_sector_databundle.log",
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_sector_databundle.py"


if config["enable"]["retrieve"] and (
    config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]
):
    datafiles = [
        "IGGIELGN_LNGs.geojson",
        "IGGIELGN_BorderPoints.geojson",
        "IGGIELGN_Productions.geojson",
        "IGGIELGN_PipeSegments.geojson",
    ]

    rule retrieve_gas_infrastructure_data:
        output:
            protected(
                expand("data/gas_network/scigrid-gas/data/{files}", files=datafiles)
            ),
        log:
            LOGS + "retrieve_gas_infrastructure_data.log",
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_gas_infrastructure_data.py"


if config["enable"]["retrieve"]:

    rule retrieve_electricity_demand:
        input:
            HTTP.remote(
                "data.open-power-system-data.org/time_series/{version}/time_series_60min_singleindex.csv".format(
                version="2019-06-05"
                    if config["snapshots"]["end"] < "2019"
                    else "2020-10-06"
                ),
                keep_local=True,
                static=True,
            ),
        output:
            RESOURCES + "load_raw.csv",
        log:
            LOGS + "retrieve_electricity_demand.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    rule retrieve_ship_raster:
        input:
            HTTP.remote(
                "https://zenodo.org/record/6953563/files/shipdensity_global.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/shipdensity_global.zip"),
        log:
            LOGS + "retrieve_ship_raster.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:

    # Downloading Copernicus Global Land Cover for land cover and land use:
    # Website: https://land.copernicus.eu/global/products/lc
    rule download_copernicus_land_cover:
        input:
            HTTP.remote(
                "zenodo.org/record/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
                static=True,
            ),
        output:
            RESOURCES
            + "Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        run:
            move(input[0], output[0])


if config["enable"]["retrieve"]:
    current_month = datetime.now().strftime("%b")
    current_year = datetime.now().strftime("%Y")
    bYYYY = f"{current_month}{current_year}"

    def check_file_exists(url):
        response = requests.head(url)
        return response.status_code == 200

    url = f"https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_{bYYYY}_Public.zip"

    if not check_file_exists(url):
        prev_month = (datetime.now() - timedelta(30)).strftime("%b")
        bYYYY = f"{prev_month}{current_year}"
        assert check_file_exists(
            f"https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_{bYYYY}_Public.zip"
        ), "The file does not exist."

    # Downloading protected area database from WDPA
    # extract the main zip and then merge the contained 3 zipped shapefiles
    # Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
    rule download_wdpa:
        input:
            HTTP.remote(
                f"d1gam3xoknrgr2.cloudfront.net/current/WDPA_{bYYYY}_Public_shp.zip",
                static=True,
                keep_local=True,
            ),
        params:
            zip=RESOURCES + f"WDPA_{bYYYY}_shp.zip",
            folder=directory(RESOURCES + f"WDPA_{bYYYY}"),
        output:
            gpkg=RESOURCES + f"WDPA_{bYYYY}.gpkg",
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
            HTTP.remote(
                f"d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{bYYYY}_Public_marine_shp.zip",
                static=True,
                keep_local=True,
            ),
        params:
            zip=RESOURCES + f"WDPA_WDOECM_{bYYYY}_marine.zip",
            folder=directory(RESOURCES + f"WDPA_WDOECM_{bYYYY}_marine"),
        output:
            gpkg=RESOURCES + f"WDPA_WDOECM_{bYYYY}_marine.gpkg",
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
            HTTP.remote(
                "https://www.eex.com/fileadmin/EEX/Downloads/EUA_Emission_Spot_Primary_Market_Auction_Report/Archive_Reports/emission-spot-primary-market-auction-report-2019-data.xls",
                keep_local=True,
                static=True,
            ),
        output:
            "data/validation/emission-spot-primary-market-auction-report-2019-data.xls",
        log:
            LOGS + "retrieve_monthly_co2_prices.log",
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
            LOGS + "retrieve_monthly_fuel_prices.log",
        resources:
            mem_mb=5000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_monthly_fuel_prices.py"
