# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

if config["enable"].get("prepare_links_p_nom", False):

    rule prepare_links_p_nom:
        output:
            "data/links_p_nom.csv",
        log:
            LOGS + "prepare_links_p_nom.log",
        threads: 1
        resources:
            mem_mb=1500,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/prepare_links_p_nom.py"


rule build_electricity_demand:
    params:
        snapshots=config["snapshots"],
        countries=config["countries"],
        load=config["load"],
    input:
        ancient("data/load_raw.csv"),
    output:
        RESOURCES + "load.csv",
    log:
        LOGS + "build_electricity_demand.log",
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_electricity_demand.py"


rule build_powerplants:
    params:
        powerplants_filter=config["electricity"]["powerplants_filter"],
        custom_powerplants=config["electricity"]["custom_powerplants"],
        countries=config["countries"],
    input:
        base_network=RESOURCES + "networks/base.nc",
        custom_powerplants="data/custom_powerplants.csv",
    output:
        RESOURCES + "powerplants.csv",
    log:
        LOGS + "build_powerplants.log",
    threads: 1
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_powerplants.py"


rule base_network:
    params:
        countries=config["countries"],
        snapshots=config["snapshots"],
        lines=config["lines"],
        links=config["links"],
        transformers=config["transformers"],
    input:
        eg_buses="data/entsoegridkit/buses.csv",
        eg_lines="data/entsoegridkit/lines.csv",
        eg_links="data/entsoegridkit/links.csv",
        eg_converters="data/entsoegridkit/converters.csv",
        eg_transformers="data/entsoegridkit/transformers.csv",
        parameter_corrections="data/parameter_corrections.yaml",
        links_p_nom="data/links_p_nom.csv",
        links_tyndp="data/links_tyndp.csv",
        country_shapes=RESOURCES + "country_shapes.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        europe_shape=RESOURCES + "europe_shape.geojson",
    output:
        RESOURCES + "networks/base.nc",
    log:
        LOGS + "base_network.log",
    benchmark:
        BENCHMARKS + "base_network"
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/base_network.py"


rule build_shapes:
    params:
        countries=config["countries"],
    input:
        naturalearth=ancient("data/bundle/naturalearth/ne_10m_admin_0_countries.shp"),
        eez=ancient("data/bundle/eez/World_EEZ_v8_2014.shp"),
        nuts3=ancient("data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp"),
        nuts3pop=ancient("data/bundle/nama_10r_3popgdp.tsv.gz"),
        nuts3gdp=ancient("data/bundle/nama_10r_3gdp.tsv.gz"),
        ch_cantons=ancient("data/bundle/ch_cantons.csv"),
        ch_popgdp=ancient("data/bundle/je-e-21.03.02.xls"),
    output:
        country_shapes=RESOURCES + "country_shapes.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        europe_shape=RESOURCES + "europe_shape.geojson",
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
    log:
        LOGS + "build_shapes.log",
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_shapes.py"


rule build_bus_regions:
    params:
        countries=config["countries"],
    input:
        country_shapes=RESOURCES + "country_shapes.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        base_network=RESOURCES + "networks/base.nc",
    output:
        regions_onshore=RESOURCES + "regions_onshore.geojson",
        regions_offshore=RESOURCES + "regions_offshore.geojson",
    log:
        LOGS + "build_bus_regions.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_bus_regions.py"


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        params:
            snapshots=config["snapshots"],
            cutouts=config["atlite"]["cutouts"],
        input:
            regions_onshore=RESOURCES + "regions_onshore.geojson",
            regions_offshore=RESOURCES + "regions_offshore.geojson",
        output:
            protected("cutouts/" + CDIR + "{cutout}.nc"),
        log:
            "logs/" + CDIR + "build_cutout/{cutout}.log",
        benchmark:
            "benchmarks/" + CDIR + "build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources:
            mem_mb=ATLITE_NPROCESSES * 1000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_cutout.py"


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        input:
            natura=ancient("data/bundle/natura/Natura2000_end2015.shp"),
            cutouts=expand("cutouts/" + CDIR + "{cutouts}.nc", **config["atlite"]),
        output:
            RESOURCES + "natura.tiff",
        resources:
            mem_mb=5000,
        log:
            LOGS + "build_natura_raster.log",
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_natura_raster.py"


rule build_ship_raster:
    input:
        ship_density="data/shipdensity_global.zip",
        cutouts=expand(
            "cutouts/" + CDIR + "{cutout}.nc",
            cutout=[
                config["renewable"][k]["cutout"]
                for k in config["electricity"]["renewable_carriers"]
            ],
        ),
    output:
        RESOURCES + "shipdensity_raster.tif",
    log:
        LOGS + "build_ship_raster.log",
    resources:
        mem_mb=5000,
    benchmark:
        BENCHMARKS + "build_ship_raster"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_ship_raster.py"


rule build_renewable_profiles:
    params:
        renewable=config["renewable"],
    input:
        base_network=RESOURCES + "networks/base.nc",
        corine=ancient("data/bundle/corine/g250_clc06_V18_5.tif"),
        natura=lambda w: (
            RESOURCES + "natura.tiff"
            if config["renewable"][w.technology]["natura"]
            else []
        ),
        gebco=ancient(
            lambda w: (
                "data/bundle/GEBCO_2014_2D.nc"
                if (
                    config["renewable"][w.technology].get("max_depth")
                    or config["renewable"][w.technology].get("min_depth")
                )
                else []
            )
        ),
        ship_density=lambda w: (
            RESOURCES + "shipdensity_raster.tif"
            if "ship_threshold" in config["renewable"][w.technology].keys()
            else []
        ),
        country_shapes=RESOURCES + "country_shapes.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        regions=lambda w: (
            RESOURCES + "regions_onshore.geojson"
            if w.technology in ("onwind", "solar")
            else RESOURCES + "regions_offshore.geojson"
        ),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config["renewable"][w.technology]["cutout"]
        + ".nc",
    output:
        profile=RESOURCES + "profile_{technology}.nc",
    log:
        LOGS + "build_renewable_profile_{technology}.log",
    benchmark:
        BENCHMARKS + "build_renewable_profiles_{technology}"
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=ATLITE_NPROCESSES * 5000,
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
        co2_price=RESOURCES + "co2_price.csv",
        fuel_price=RESOURCES + "monthly_fuel_price.csv",
    log:
        LOGS + "build_monthly_prices.log",
    threads: 1
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_monthly_prices.py"


rule build_hydro_profile:
    params:
        hydro=config["renewable"]["hydro"],
        countries=config["countries"],
    input:
        country_shapes=RESOURCES + "country_shapes.geojson",
        eia_hydro_generation="data/eia_hydro_annual_generation.csv",
        cutout=f"cutouts/" + CDIR + config["renewable"]["hydro"]["cutout"] + ".nc",
    output:
        RESOURCES + "profile_hydro.nc",
    log:
        LOGS + "build_hydro_profile.log",
    resources:
        mem_mb=5000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_hydro_profile.py"


if config["lines"]["dynamic_line_rating"]["activate"]:

    rule build_line_rating:
        input:
            base_network=RESOURCES + "networks/base.nc",
            cutout="cutouts/"
            + CDIR
            + config["lines"]["dynamic_line_rating"]["cutout"]
            + ".nc",
        output:
            output=RESOURCES + "networks/line_rating.nc",
        log:
            LOGS + "build_line_rating.log",
        benchmark:
            BENCHMARKS + "build_line_rating"
        threads: ATLITE_NPROCESSES
        resources:
            mem_mb=ATLITE_NPROCESSES * 1000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_line_rating.py"


rule add_electricity:
    params:
        length_factor=config["lines"]["length_factor"],
        scaling_factor=config["load"]["scaling_factor"],
        countries=config["countries"],
        renewable=config["renewable"],
        electricity=config["electricity"],
        conventional=config["conventional"],
        costs=config["costs"],
    input:
        **{
            f"profile_{tech}": RESOURCES + f"profile_{tech}.nc"
            for tech in config["electricity"]["renewable_carriers"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            if carrier in config["electricity"]["conventional_carriers"]
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network=RESOURCES + "networks/base.nc",
        line_rating=RESOURCES + "networks/line_rating.nc"
        if config["lines"]["dynamic_line_rating"]["activate"]
        else RESOURCES + "networks/base.nc",
        tech_costs=COSTS,
        regions=RESOURCES + "regions_onshore.geojson",
        powerplants=RESOURCES + "powerplants.csv",
        hydro_capacities=ancient("data/bundle/hydro_capacities.csv"),
        geth_hydro_capacities="data/geth2015_hydro_capacities.csv",
        unit_commitment="data/unit_commitment.csv",
        fuel_price=RESOURCES + "monthly_fuel_price.csv"
        if config["conventional"]["dynamic_fuel_price"]
        else [],
        load=RESOURCES + "load.csv",
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
    output:
        RESOURCES + "networks/elec.nc",
    log:
        LOGS + "add_electricity.log",
    benchmark:
        BENCHMARKS + "add_electricity"
    threads: 1
    resources:
        mem_mb=10000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_electricity.py"


rule simplify_network:
    params:
        simplify_network=config["clustering"]["simplify_network"],
        aggregation_strategies=config["clustering"].get("aggregation_strategies", {}),
        focus_weights=config.get("focus_weights", None),
        renewable_carriers=config["electricity"]["renewable_carriers"],
        max_hours=config["electricity"]["max_hours"],
        length_factor=config["lines"]["length_factor"],
        p_max_pu=config["links"].get("p_max_pu", 1.0),
        costs=config["costs"],
    input:
        network=RESOURCES + "networks/elec.nc",
        tech_costs=COSTS,
        regions_onshore=RESOURCES + "regions_onshore.geojson",
        regions_offshore=RESOURCES + "regions_offshore.geojson",
    output:
        network=RESOURCES + "networks/elec_s{simpl}.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}.geojson",
        regions_offshore=RESOURCES + "regions_offshore_elec_s{simpl}.geojson",
        busmap=RESOURCES + "busmap_elec_s{simpl}.csv",
        connection_costs=RESOURCES + "connection_costs_s{simpl}.csv",
    log:
        LOGS + "simplify_network/elec_s{simpl}.log",
    benchmark:
        BENCHMARKS + "simplify_network/elec_s{simpl}"
    threads: 1
    resources:
        mem_mb=12000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/simplify_network.py"


rule cluster_network:
    params:
        cluster_network=config["clustering"]["cluster_network"],
        aggregation_strategies=config["clustering"].get("aggregation_strategies", {}),
        custom_busmap=config["enable"].get("custom_busmap", False),
        focus_weights=config.get("focus_weights", None),
        renewable_carriers=config["electricity"]["renewable_carriers"],
        conventional_carriers=config["electricity"].get("conventional_carriers", []),
        max_hours=config["electricity"]["max_hours"],
        length_factor=config["lines"]["length_factor"],
        costs=config["costs"],
    input:
        network=RESOURCES + "networks/elec_s{simpl}.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}.geojson",
        regions_offshore=RESOURCES + "regions_offshore_elec_s{simpl}.geojson",
        busmap=ancient(RESOURCES + "busmap_elec_s{simpl}.csv"),
        custom_busmap=(
            "data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            if config["enable"].get("custom_busmap", False)
            else []
        ),
        tech_costs=COSTS,
    output:
        network=RESOURCES + "networks/elec_s{simpl}_{clusters}.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore=RESOURCES + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
        busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
        linemap=RESOURCES + "linemap_elec_s{simpl}_{clusters}.csv",
    log:
        LOGS + "cluster_network/elec_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "cluster_network/elec_s{simpl}_{clusters}"
    threads: 1
    resources:
        mem_mb=10000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/cluster_network.py"


rule add_extra_components:
    params:
        extendable_carriers=config["electricity"]["extendable_carriers"],
        max_hours=config["electricity"]["max_hours"],
        costs=config["costs"],
    input:
        network=RESOURCES + "networks/elec_s{simpl}_{clusters}.nc",
        tech_costs=COSTS,
    output:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc",
    log:
        LOGS + "add_extra_components/elec_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "add_extra_components/elec_s{simpl}_{clusters}_ec"
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_extra_components.py"


rule prepare_network:
    params:
        links=config["links"],
        lines=config["lines"],
        co2base=config["electricity"]["co2base"],
        co2limit=config["electricity"]["co2limit"],
        gaslimit=config["electricity"].get("gaslimit"),
        max_hours=config["electricity"]["max_hours"],
        costs=config["costs"],
    input:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc",
        tech_costs=COSTS,
        co2_price=lambda w: RESOURCES + "co2_price.csv" if "Ept" in w.opts else [],
    output:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        LOGS + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log",
    benchmark:
        (BENCHMARKS + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}")
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/prepare_network.py"
