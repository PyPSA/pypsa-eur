# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

if config["enable"].get("retrieve_databundle", True):
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
            expand("data/bundle/{file}", file=datafiles),
        log:
            LOGS + "retrieve_databundle.log",
        resources:
            mem_mb=1000,
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_databundle.py"


if config["enable"].get("retrieve_cutout", True):

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


if config["enable"].get("retrieve_cost_data", True):

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


if config["enable"].get("retrieve_natura_raster", True):

    rule retrieve_natura_raster:
        input:
            HTTP.remote(
                "zenodo.org/record/4706686/files/natura.tiff",
                keep_local=True,
                static=True,
            ),
        output:
            RESOURCES + "natura.tiff",
        log:
            LOGS + "retrieve_natura_raster.log",
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


if config["enable"].get("retrieve_sector_databundle", True):
    datafiles = [
        "data/eea/UNFCCC_v23.csv",
        "data/switzerland-sfoe/switzerland-new_format.csv",
        "data/nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson",
        "data/myb1-2017-nitro.xls",
        "data/Industrial_Database.csv",
        "data/emobility/KFZ__count",
        "data/emobility/Pkw__count",
        "data/h2_salt_caverns_GWh_per_sqkm.geojson",
        directory("data/eurostat-energy_balances-june_2016_edition"),
        directory("data/eurostat-energy_balances-may_2018_edition"),
        directory("data/jrc-idees-2015"),
    ]

    rule retrieve_sector_databundle:
        output:
            *datafiles,
        log:
            LOGS + "retrieve_sector_databundle.log",
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_sector_databundle.py"


if config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]:
    datafiles = [
        "IGGIELGN_LNGs.geojson",
        "IGGIELGN_BorderPoints.geojson",
        "IGGIELGN_Productions.geojson",
        "IGGIELGN_PipeSegments.geojson",
    ]

    rule retrieve_gas_infrastructure_data:
        output:
            expand("data/gas_network/scigrid-gas/data/{files}", files=datafiles),
        log:
            LOGS + "retrieve_gas_infrastructure_data.log",
        retries: 2
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/retrieve_gas_infrastructure_data.py"


rule retrieve_electricity_demand:
    input:
        HTTP.remote(
            "data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv",
            keep_local=True,
            static=True,
        ),
    output:
        "data/load_raw.csv",
    log:
        LOGS + "retrieve_electricity_demand.log",
    resources:
        mem_mb=5000,
    retries: 2
    run:
        move(input[0], output[0])


rule retrieve_ship_raster:
    input:
        HTTP.remote(
            "https://zenodo.org/record/6953563/files/shipdensity_global.zip",
            keep_local=True,
            static=True,
        ),
    output:
        "data/shipdensity_global.zip",
    log:
        LOGS + "retrieve_ship_raster.log",
    resources:
        mem_mb=5000,
    retries: 2
    run:
        move(input[0], output[0])
