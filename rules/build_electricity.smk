# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

if config["enable"].get("prepare_links_p_nom", False):

    rule prepare_links_p_nom:
        output:
            "data/links_p_nom.csv",
        log:
            logs("prepare_links_p_nom.log"),
        threads: 1
        resources:
            mem_mb=1500,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/prepare_links_p_nom.py"


rule build_electricity_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        countries=config_provider("countries"),
        load=config_provider("load"),
    input:
        reported=ancient("data/electricity_demand_raw.csv"),
        synthetic=lambda w: (
            ancient("data/load_synthetic_raw.csv")
            if config_provider("load", "supplement_synthetic")(w)
            else []
        ),
    output:
        resources("electricity_demand.csv"),
    log:
        logs("build_electricity_demand.log"),
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_electricity_demand.py"


rule build_powerplants:
    params:
        powerplants_filter=config_provider("electricity", "powerplants_filter"),
        custom_powerplants=config_provider("electricity", "custom_powerplants"),
        everywhere_powerplants=config_provider("electricity", "everywhere_powerplants"),
        countries=config_provider("countries"),
    input:
        base_network=resources("networks/base.nc"),
        custom_powerplants="data/custom_powerplants.csv",
    output:
        resources("powerplants.csv"),
    log:
        logs("build_powerplants.log"),
    threads: 1
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_powerplants.py"


rule base_network:
    params:
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        line_projects=config_provider("line_projects"),
        lines=config_provider("lines"),
        transformers=config_provider("transformers"),
    input:
        eg_buses="data/entsoegridkit/buses.csv",
        eg_lines="data/entsoegridkit/lines.csv",
        eg_links="data/entsoegridkit/links.csv",
        eg_converters="data/entsoegridkit/converters.csv",
        eg_transformers="data/entsoegridkit/transformers.csv",
        parameter_corrections="data/parameter_corrections.yaml",
        links_p_nom="data/links_p_nom.csv",
        links_tyndp="data/links_tyndp.csv",
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        europe_shape=resources("europe_shape.geojson"),
    output:
        base_network=resources("networks/base.nc"),
        regions_onshore=resources("regions_onshore.geojson"),
        regions_offshore=resources("regions_offshore.geojson"),
    log:
        logs("base_network.log"),
    benchmark:
        benchmarks("base_network")
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/base_network.py"


rule build_shapes:
    params:
        countries=config_provider("countries"),
    input:
        naturalearth=ancient("data/bundle/naturalearth/ne_10m_admin_0_countries.shp"),
        eez=ancient("data/bundle/eez/World_EEZ_v8_2014.shp"),
        nuts3=ancient("data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp"),
        nuts3pop=ancient("data/bundle/nama_10r_3popgdp.tsv.gz"),
        nuts3gdp=ancient("data/bundle/nama_10r_3gdp.tsv.gz"),
        ch_cantons=ancient("data/ch_cantons.csv"),
        ch_popgdp=ancient("data/bundle/je-e-21.03.02.xls"),
    output:
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        europe_shape=resources("europe_shape.geojson"),
        nuts3_shapes=resources("nuts3_shapes.geojson"),
    log:
        logs("build_shapes.log"),
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_shapes.py"


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        params:
            snapshots=config_provider("snapshots"),
            cutouts=config_provider("atlite", "cutouts"),
        input:
            regions_onshore=resources("regions_onshore.geojson"),
            regions_offshore=resources("regions_offshore.geojson"),
        output:
            protected("cutouts/" + CDIR + "{cutout}.nc"),
        log:
            logs(CDIR + "build_cutout/{cutout}.log"),
        benchmark:
            "benchmarks/" + CDIR + "build_cutout_{cutout}"
        threads: config["atlite"].get("nprocesses", 4)
        resources:
            mem_mb=config["atlite"].get("nprocesses", 4) * 1000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_cutout.py"


rule build_ship_raster:
    input:
        ship_density="data/shipdensity_global.zip",
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("atlite", "default_cutout")(w)
        + ".nc",
    output:
        resources("shipdensity_raster.tif"),
    log:
        logs("build_ship_raster.log"),
    resources:
        mem_mb=5000,
    benchmark:
        benchmarks("build_ship_raster")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_ship_raster.py"


rule determine_availability_matrix_MD_UA:
    input:
        copernicus="data/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        wdpa="data/WDPA.gpkg",
        wdpa_marine="data/WDPA_WDOECM_marine.gpkg",
        gebco=lambda w: (
            "data/bundle/gebco/GEBCO_2014_2D.nc"
            if config_provider("renewable", w.technology)(w).get("max_depth")
            else []
        ),
        ship_density=lambda w: (
            resources("shipdensity_raster.tif")
            if "ship_threshold" in config_provider("renewable", w.technology)(w).keys()
            else []
        ),
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        regions=lambda w: (
            resources("regions_onshore.geojson")
            if w.technology in ("onwind", "solar")
            else resources("regions_offshore.geojson")
        ),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("renewable", w.technology, "cutout")(w)
        + ".nc",
    output:
        availability_matrix=resources("availability_matrix_MD-UA_{technology}.nc"),
        availability_map=resources("availability_matrix_MD-UA_{technology}.png"),
    log:
        logs("determine_availability_matrix_MD_UA_{technology}.log"),
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/determine_availability_matrix_MD_UA.py"


# Optional input when having Ukraine (UA) or Moldova (MD) in the countries list
def input_ua_md_availability_matrix(w):
    countries = set(config_provider("countries")(w))
    if {"UA", "MD"}.intersection(countries):
        return {
            "availability_matrix_MD_UA": resources(
                "availability_matrix_MD-UA_{technology}.nc"
            )
        }
    return {}


rule build_renewable_profiles:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        renewable=config_provider("renewable"),
    input:
        unpack(input_ua_md_availability_matrix),
        base_network=resources("networks/base.nc"),
        corine=ancient("data/bundle/corine/g250_clc06_V18_5.tif"),
        natura=lambda w: (
            "data/bundle/natura/natura.tiff"
            if config_provider("renewable", w.technology, "natura")(w)
            else []
        ),
        luisa=lambda w: (
            "data/LUISA_basemap_020321_50m.tif"
            if config_provider("renewable", w.technology, "luisa")(w)
            else []
        ),
        gebco=ancient(
            lambda w: (
                "data/bundle/gebco/GEBCO_2014_2D.nc"
                if (
                    config_provider("renewable", w.technology)(w).get("max_depth")
                    or config_provider("renewable", w.technology)(w).get("min_depth")
                )
                else []
            )
        ),
        ship_density=lambda w: (
            resources("shipdensity_raster.tif")
            if "ship_threshold" in config_provider("renewable", w.technology)(w).keys()
            else []
        ),
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        regions=lambda w: (
            resources("regions_onshore.geojson")
            if w.technology in ("onwind", "solar")
            else resources("regions_offshore.geojson")
        ),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("renewable", w.technology, "cutout")(w)
        + ".nc",
    output:
        profile=resources("profile_{technology}.nc"),
    log:
        logs("build_renewable_profile_{technology}.log"),
    benchmark:
        benchmarks("build_renewable_profiles_{technology}")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 5000,
    wildcard_constraints:
        technology="(?!hydro).*",  # Any technology other than hydro
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_renewable_profiles.py"


rule build_monthly_prices:
    input:
        co2_price_raw="data/validation/emission-spot-primary-market-auction-report-2019-data.xls",
        fuel_price_raw="data/validation/energy-price-trends-xlsx-5619002.xlsx",
    output:
        co2_price=resources("co2_price.csv"),
        fuel_price=resources("monthly_fuel_price.csv"),
    log:
        logs("build_monthly_prices.log"),
    threads: 1
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_monthly_prices.py"


rule build_hydro_profile:
    params:
        hydro=config_provider("renewable", "hydro"),
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        country_shapes=resources("country_shapes.geojson"),
        eia_hydro_generation="data/eia_hydro_annual_generation.csv",
        eia_hydro_capacity="data/eia_hydro_annual_capacity.csv",
        era5_runoff="data/era5-annual-runoff-per-country.csv",
        cutout=lambda w: f"cutouts/"
        + CDIR
        + config_provider("renewable", "hydro", "cutout")(w)
        + ".nc",
    output:
        profile=resources("profile_hydro.nc"),
    log:
        logs("build_hydro_profile.log"),
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_hydro_profile.py"


rule build_line_rating:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        base_network=resources("networks/base.nc"),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("lines", "dynamic_line_rating", "cutout")(w)
        + ".nc",
    output:
        output=resources("networks/line_rating.nc"),
    log:
        logs("build_line_rating.log"),
    benchmark:
        benchmarks("build_line_rating")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 1000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_line_rating.py"


def input_profile_tech(w):
    return {
        f"profile_{tech}": resources(f"profile_{tech}.nc")
        for tech in config_provider("electricity", "renewable_carriers")(w)
    }


def input_conventional(w):
    return {
        f"conventional_{carrier}_{attr}": fn
        for carrier, d in config_provider("conventional", default={None: {}})(w).items()
        if carrier in config_provider("electricity", "conventional_carriers")(w)
        for attr, fn in d.items()
        if str(fn).startswith("data/")
    }


rule add_electricity:
    params:
        length_factor=config_provider("lines", "length_factor"),
        scaling_factor=config_provider("load", "scaling_factor"),
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        renewable=config_provider("renewable"),
        electricity=config_provider("electricity"),
        conventional=config_provider("conventional"),
        costs=config_provider("costs"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        unpack(input_profile_tech),
        unpack(input_conventional),
        base_network=resources("networks/base.nc"),
        line_rating=lambda w: (
            resources("networks/line_rating.nc")
            if config_provider("lines", "dynamic_line_rating", "activate")(w)
            else resources("networks/base.nc")
        ),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
        regions=resources("regions_onshore.geojson"),
        powerplants=resources("powerplants.csv"),
        hydro_capacities=ancient("data/hydro_capacities.csv"),
        geth_hydro_capacities="data/geth2015_hydro_capacities.csv",
        unit_commitment="data/unit_commitment.csv",
        fuel_price=lambda w: (
            resources("monthly_fuel_price.csv")
            if config_provider("conventional", "dynamic_fuel_price")(w)
            else []
        ),
        load=resources("electricity_demand.csv"),
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        ua_md_gdp="data/GDP_PPP_30arcsec_v3_mapped_default.csv",
    output:
        resources("networks/elec.nc"),
    log:
        logs("add_electricity.log"),
    benchmark:
        benchmarks("add_electricity")
    threads: 1
    resources:
        mem_mb=10000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_electricity.py"


rule simplify_network:
    params:
        simplify_network=config_provider("clustering", "simplify_network"),
        aggregation_strategies=config_provider(
            "clustering", "aggregation_strategies", default={}
        ),
        focus_weights=config_provider("clustering", "focus_weights", default=None),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        max_hours=config_provider("electricity", "max_hours"),
        length_factor=config_provider("lines", "length_factor"),
        p_max_pu=config_provider("links", "p_max_pu", default=1.0),
        costs=config_provider("costs"),
    input:
        network=resources("networks/elec.nc"),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
        regions_onshore=resources("regions_onshore.geojson"),
        regions_offshore=resources("regions_offshore.geojson"),
    output:
        network=resources("networks/elec_s{simpl}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}.geojson"),
        busmap=resources("busmap_elec_s{simpl}.csv"),
    log:
        logs("simplify_network/elec_s{simpl}.log"),
    benchmark:
        benchmarks("simplify_network/elec_s{simpl}")
    threads: 1
    resources:
        mem_mb=12000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/simplify_network.py"


rule cluster_network:
    params:
        cluster_network=config_provider("clustering", "cluster_network"),
        aggregation_strategies=config_provider(
            "clustering", "aggregation_strategies", default={}
        ),
        custom_busmap=config_provider("enable", "custom_busmap", default=False),
        focus_weights=config_provider("clustering", "focus_weights", default=None),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        conventional_carriers=config_provider(
            "electricity", "conventional_carriers", default=[]
        ),
        max_hours=config_provider("electricity", "max_hours"),
        length_factor=config_provider("lines", "length_factor"),
        costs=config_provider("costs"),
    input:
        network=resources("networks/elec_s{simpl}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}.geojson"),
        busmap=ancient(resources("busmap_elec_s{simpl}.csv")),
        custom_busmap=lambda w: (
            "data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            if config_provider("enable", "custom_busmap", default=False)(w)
            else []
        ),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
    output:
        network=resources("networks/elec_s{simpl}_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
        busmap=resources("busmap_elec_s{simpl}_{clusters}.csv"),
        linemap=resources("linemap_elec_s{simpl}_{clusters}.csv"),
    log:
        logs("cluster_network/elec_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("cluster_network/elec_s{simpl}_{clusters}")
    threads: 1
    resources:
        mem_mb=10000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/cluster_network.py"


rule add_extra_components:
    params:
        extendable_carriers=config_provider("electricity", "extendable_carriers"),
        max_hours=config_provider("electricity", "max_hours"),
        costs=config_provider("costs"),
    input:
        network=resources("networks/elec_s{simpl}_{clusters}.nc"),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
    output:
        resources("networks/elec_s{simpl}_{clusters}_ec.nc"),
    log:
        logs("add_extra_components/elec_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("add_extra_components/elec_s{simpl}_{clusters}_ec")
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_extra_components.py"


rule prepare_network:
    params:
        time_resolution=config_provider("clustering", "temporal", "resolution_elec"),
        links=config_provider("links"),
        lines=config_provider("lines"),
        co2base=config_provider("electricity", "co2base"),
        co2limit_enable=config_provider("electricity", "co2limit_enable", default=False),
        co2limit=config_provider("electricity", "co2limit"),
        gaslimit_enable=config_provider("electricity", "gaslimit_enable", default=False),
        gaslimit=config_provider("electricity", "gaslimit"),
        max_hours=config_provider("electricity", "max_hours"),
        costs=config_provider("costs"),
        adjustments=config_provider("adjustments", "electricity"),
        autarky=config_provider("electricity", "autarky", default={}),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        resources("networks/elec_s{simpl}_{clusters}_ec.nc"),
        tech_costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}.csv"
        ),
        co2_price=lambda w: resources("co2_price.csv") if "Ept" in w.opts else [],
    output:
        resources("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
    log:
        logs("prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log"),
    benchmark:
        (benchmarks("prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"))
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/prepare_network.py"
