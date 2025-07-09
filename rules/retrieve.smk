# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import requests
from datetime import datetime, timedelta
from shutil import move, unpack_archive
from shutil import copy as shcopy
from zipfile import ZipFile

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
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_databundle.py"


if (EUROSTAT_BALANCES_DATASET := dataset_version("eurostat_balances"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eurostat_balances:
        params:
            url=EUROSTAT_BALANCES_DATASET["url"],
        output:
            zip=f"{EUROSTAT_BALANCES_DATASET['folder']}/balances.zip",
            directory=directory(f"{EUROSTAT_BALANCES_DATASET['folder']}"),
        run:
            import os
            import requests
            from zipfile import ZipFile
            from pathlib import Path

            response = requests.get(params["url"])
            with open(output.zip, "wb") as f:
                f.write(response.content)

            output_folder = Path(output["zip"]).parent
            unpack_archive(output.zip, output_folder)



if config["enable"]["retrieve"] and config["enable"].get("retrieve_jrc_idees", True):

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
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_eurostat_household_data.py"


if (EU_NUTS2013_DATASET := dataset_version("eu_nuts2013"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_eu_nuts_2013:
        input:
            shapes=storage(EU_NUTS2013_DATASET["url"]),
        output:
            shapes_level_3=f"{EU_NUTS2013_DATASET["folder"]}/NUTS_RG_03M_2013_4326_LEVL_3.geojson",
            shapes_level_2=f"{EU_NUTS2013_DATASET["folder"]}/NUTS_RG_03M_2013_4326_LEVL_2.geojson",
        params:
            zip_file=f"{EU_NUTS2013_DATASET['folder']}/ref-nuts-2013-03m.geojson.zip",
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
            shapes_level_3=f"{EU_NUTS2021_DATASET["folder"]}/NUTS_RG_01M_2021_4326_LEVL_3.geojson",
            shapes_level_2=f"{EU_NUTS2021_DATASET["folder"]}/NUTS_RG_01M_2021_4326_LEVL_2.geojson",
            shapes_level_1=f"{EU_NUTS2021_DATASET["folder"]}/NUTS_RG_01M_2021_4326_LEVL_1.geojson",
            shapes_level_0=f"{EU_NUTS2021_DATASET["folder"]}/NUTS_RG_01M_2021_4326_LEVL_0.geojson",
        params:
            zip_file=f"{EU_NUTS2021_DATASET['folder']}/ref-nuts-2021-01m.geojson.zip",
        run:
            os.rename(input.shapes, params.zip_file)
            with ZipFile(params.zip_file, "r") as zip_ref:
                for level in ["LEVL_3", "LEVL_2", "LEVL_1", "LEVL_0"]:
                    filename = f"NUTS_RG_01M_2021_4326_{level}.geojson"
                    zip_ref.extract(filename, Path(output.shapes_level_0).parent)
                    extracted_file = Path(output.shapes_level_0).parent / filename
                    extracted_file.rename(
                        getattr(output, f"shapes_level_{level[-1]}")
                    )



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


if config["enable"]["retrieve"] and config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            storage(
                "https://zenodo.org/records/14936211/files/{cutout}.nc",
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
            "../envs/environment.yaml"
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
            "../envs/environment.yaml"
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
            "../envs/environment.yaml"
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


if (JRC_ENSPRESO_BIOMASS_DATASET := dataset_version("jrc_enspreso_biomass"))[
    "source"
] in [
    "primary",
    "archive",
]:

    rule retrieve_jrc_enspreso_biomass:
        input:
            storage(
                f"{JRC_ENSPRESO_BIOMASS_DATASET["url"]}",
                keep_local=True,
            ),
        output:
            f"{JRC_ENSPRESO_BIOMASS_DATASET["folder"]}/ENSPRESO_BIOMASS.xlsx",
        retries: 1
        run:
            move(input[0], output[0])


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
                keep_local=True,
            ),
        output:
            f"{HOTMAPS_INDUSTRIAL_SITES["folder"]}/Industrial_Database.csv",
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
            f"{NITROGEN_STATISTICS_DATASET['folder']}/nitro-ert.xlsx",
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
            import os
            import requests

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



if (CO2STOP_DATASET := dataset_version("co2stop"))["source"] in [
    "primary",
    "archive",
]:

    rule retrieve_co2stop:
        output:
            storage_table=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Storage_Units.csv",
            storage_map=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_Polygons Data/StorageUnits_March13.kml",
            traps_table1=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps.csv",
            traps_table2=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps_Temp.csv",
            traps_table3=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_DataInterrogationSystem/Hydrocarbon_Traps1.csv",
            traps_map=f"{CO2STOP_DATASET['folder']}/CO2JRC_OpenFormats/CO2Stop_Polygons Data/DaughterUnits_March13.kml",
        run:
            import requests

            response = requests.get(
                CO2STOP_DATASET["url"],
            )
            zip_file = f"{CO2STOP_DATASET['folder']}/co2jrc_openformats.zip"
            output_folder = CO2STOP_DATASET["folder"]

            with open(zip_file, "wb") as f:
                f.write(response.content)
            unpack_archive(zip_file, output_folder)



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
            os.rename(input.xlsx, output.xlsx)


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
            zip=storage(url, keep_local=True),
        params:
            zip="data/WDPA_shp.zip",
            folder=directory("data/WDPA"),
        output:
            gpkg="data/WDPA.gpkg",
        run:
            shcopy(input.zip, params.zip)
            unpack_archive(params.zip, params.folder)

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
            zip=storage(
                f"https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{bYYYY}_Public_marine_shp.zip",
                keep_local=True,
            ),
        params:
            zip="data/WDPA_WDOECM_marine.zip",
            folder=directory("data/WDPA_WDOECM_marine"),
        output:
            gpkg="data/WDPA_WDOECM_marine.gpkg",
        run:
            shcopy(input.zip, params.zip)
            unpack_archive(params.zip, params.folder)

            for i in range(3):
                # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
                layer_path = f"/vsizip/{params.folder}/WDPA_WDOECM_{bYYYY}_Public_marine_shp_{i}.zip"
                print(f"Adding layer {i+1} of 3 to combined output file.")
                shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")



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
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_monthly_fuel_prices.py"


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

    OSM_URL = dataset_version("osm")["url"]

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


if config["enable"]["retrieve"]:

    rule retrieve_osm_boundaries:
        output:
            json="data/osm-boundaries/json/{country}_adm1.json",
        log:
            "logs/retrieve_osm_boundaries_{country}_adm1.log",
        threads: 1
        conda:
            "../envs/environment.yaml"
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


if config["enable"]["retrieve"]:

    rule retrieve_jrc_ardeco:
        output:
            ardeco_gdp="data/jrc-ardeco/ARDECO-SUVGDP.2021.table.csv",
            ardeco_pop="data/jrc-ardeco/ARDECO-SNPTD.2021.table.csv",
        run:
            import requests

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



if config["enable"]["retrieve"]:

    rule retrieve_aquifer_data_bgr:
        input:
            zip=storage(
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
            with ZipFile(input.zip, "r") as zip_ref:
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
