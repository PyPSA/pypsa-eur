# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule build_electricity_demand:
    message:
        "Building electricity demand time series"
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        countries=config_provider("countries"),
        load=config_provider("load"),
    input:
        opsd=rules.retrieve_electricity_demand_opsd.output["csv"],
        neso=rules.retrieve_electricity_demand_neso.output["csv"],
        entsoe=rules.retrieve_electricity_demand_entsoe.output["csv"],
        synthetic=lambda w: (
            ancient(rules.retrieve_synthetic_electricity_demand.output["csv"])
            if config_provider("load", "supplement_synthetic")(w)
            else []
        ),
    output:
        resources("electricity_demand.csv"),
    log:
        logs("build_electricity_demand.log"),
    benchmark:
        benchmarks("build_electricity_demand")
    resources:
        mem_mb=5000,
    script:
        scripts("build_electricity_demand.py")


rule build_powerplants:
    message:
        "Building powerplant list for {wildcards.clusters} clusters"
    params:
        powerplants_filter=config_provider("electricity", "powerplants_filter"),
        custom_powerplants=config_provider("electricity", "custom_powerplants"),
        everywhere_powerplants=config_provider("electricity", "everywhere_powerplants"),
        countries=config_provider("countries"),
    input:
        network=resources("networks/base_s_{clusters}.nc"),
        powerplants=rules.retrieve_powerplants.output["powerplants"],
        custom_powerplants="data/custom_powerplants.csv",
    output:
        resources("powerplants_s_{clusters}.csv"),
    log:
        logs("build_powerplants_s_{clusters}.log"),
    benchmark:
        benchmarks("build_powerplants_s_{clusters}")
    threads: 1
    resources:
        mem_mb=7000,
    script:
        scripts("build_powerplants.py")


def input_base_network(w):
    base_network = config_provider("electricity", "base_network")(w)
    source = config_provider("data", "osm", "source")(w)
    components = {"buses", "lines", "links", "converters", "transformers"}
    if (base_network == "osm") and (source == "archive"):
        OSM_DATASET = dataset_version("osm")
        inputs = {c: f"{OSM_DATASET['folder']}/{c}.csv" for c in components}
    elif base_network == "osm" and (source == "build"):
        inputs = {c: resources(f"osm/build/{c}.csv") for c in components}
    elif base_network == "tyndp":
        inputs = {c: resources(f"tyndp/build/{c}.csv") for c in components}
    elif base_network == "entsoegridkit":
        inputs = {c: f"data/entsoegridkit/{c}.csv" for c in components}
        inputs["parameter_corrections"] = "data/parameter_corrections.yaml"
        inputs["links_p_nom"] = "data/links_p_nom.csv"
    return inputs


rule base_network:
    message:
        "Building base network"
    params:
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        lines=config_provider("lines"),
        links=config_provider("links"),
        transformers=config_provider("transformers"),
        clustering=config_provider("clustering", "mode"),
        admin_levels=config_provider("clustering", "administrative"),
    input:
        unpack(input_base_network),
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        europe_shape=resources("europe_shape.geojson"),
    output:
        base_network=resources("networks/base.nc"),
        regions_onshore=resources("regions_onshore.geojson"),
        regions_offshore=resources("regions_offshore.geojson"),
        admin_shapes=resources("admin_shapes.geojson"),
    log:
        logs("base_network.log"),
    benchmark:
        benchmarks("base_network")
    threads: 4
    resources:
        mem_mb=2000,
    script:
        scripts("base_network.py")


rule build_osm_boundaries:
    message:
        "Building OSM boundaries for {wildcards.country}"
    input:
        json=f"{OSM_BOUNDARIES_DATASET['folder']}/{{country}}_adm1.json",
        eez=ancient(rules.retrieve_eez.output["gpkg"]),
    output:
        boundary=f"data/osm_boundaries/build/{OSM_BOUNDARIES_DATASET['version']}/{{country}}_adm1.geojson",
    log:
        "logs/build_osm_boundaries_{country}.log",
    threads: 1
    resources:
        mem_mb=1500,
    script:
        scripts("build_osm_boundaries.py")


rule build_bidding_zones:
    message:
        "Building bidding zones"
    params:
        countries=config_provider("countries"),
        remove_islands=config_provider(
            "clustering", "build_bidding_zones", "remove_islands"
        ),
        aggregate_to_tyndp=config_provider(
            "clustering", "build_bidding_zones", "aggregate_to_tyndp"
        ),
    input:
        bidding_zones_entsoepy=f"{BIDDING_ZONES_ENTSOEPY_DATASET['folder']}/bidding_zones_entsoepy.geojson",
        bidding_zones_electricitymaps=f"{BIDDING_ZONES_ELECTRICITYMAPS_DATASET['folder']}/bidding_zones_electricitymaps.geojson",
    output:
        file=resources("bidding_zones.geojson"),
    log:
        logs("build_bidding_zones.log"),
    threads: 1
    resources:
        mem_mb=1500,
    script:
        scripts("build_bidding_zones.py")


rule build_shapes:
    message:
        "Building geographical shapes"
    params:
        config_provider("clustering", "mode"),
        countries=config_provider("countries"),
    input:
        eez=ancient(rules.retrieve_eez.output["gpkg"]),
        nuts3_2021=rules.retrieve_eu_nuts_2021.output["shapes_level_3"],
        ba_adm1=f"data/osm_boundaries/build/{OSM_BOUNDARIES_DATASET['version']}/BA_adm1.geojson",
        md_adm1=f"data/osm_boundaries/build/{OSM_BOUNDARIES_DATASET['version']}/MD_adm1.geojson",
        ua_adm1=f"data/osm_boundaries/build/{OSM_BOUNDARIES_DATASET['version']}/UA_adm1.geojson",
        xk_adm1=f"data/osm_boundaries/build/{OSM_BOUNDARIES_DATASET['version']}/XK_adm1.geojson",
        nuts3_gdp=rules.retrieve_jrc_ardeco.output["ardeco_gdp"],
        nuts3_pop=rules.retrieve_jrc_ardeco.output["ardeco_pop"],
        bidding_zones=lambda w: (
            resources("bidding_zones.geojson")
            if config_provider("clustering", "mode")(w) == "administrative"
            else []
        ),
        other_gdp=rules.retrieve_gdp_per_capita.output["gdp"],
        other_pop=rules.retrieve_population_count.output["tif"],
    output:
        country_shapes=resources("country_shapes.geojson"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        europe_shape=resources("europe_shape.geojson"),
        nuts3_shapes=resources("nuts3_shapes.geojson"),
    log:
        logs("build_shapes.log"),
    benchmark:
        benchmarks("build_shapes")
    threads: 1
    resources:
        mem_mb=1500,
    script:
        scripts("build_shapes.py")


if CUTOUT_DATASET["source"] in ["build"]:

    rule build_cutout:
        message:
            "Building cutout data for {wildcards.cutout}"
        params:
            cutouts=config_provider("atlite", "cutouts"),
        output:
            cutout=CUTOUT_DATASET["folder"] / "{cutout}.nc",
        log:
            "logs/build_cutout/{cutout}.log",
        benchmark:
            "benchmarks/build_cutout/{cutout}"
        threads: config["atlite"].get("nprocesses", 4)
        resources:
            mem_mb=config["atlite"].get("nprocesses", 4) * 1000,
        script:
            scripts("build_cutout.py")


rule build_ship_raster:
    message:
        "Building ship density raster"
    input:
        ship_density=rules.retrieve_ship_raster.output["zip_file"],
        cutout=lambda w: input_cutout(w),
    output:
        resources("shipdensity_raster.tif"),
    log:
        logs("build_ship_raster.log"),
    resources:
        mem_mb=5000,
    benchmark:
        benchmarks("build_ship_raster")
    script:
        scripts("build_ship_raster.py")


rule determine_availability_matrix_MD_UA:
    message:
        "Determining availability matrix for {wildcards.clusters} clusters and {wildcards.technology} technology"
    params:
        renewable=config_provider("renewable"),
    input:
        copernicus=rules.download_copernicus_land_cover.output["tif"],
        wdpa=rules.retrieve_wdpa.output["gpkg"],
        wdpa_marine=rules.retrieve_wdpa_marine.output["gpkg"],
        gebco=lambda w: (
            rules.retrieve_gebco.output["gebco"]
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
            resources("regions_onshore_base_s_{clusters}.geojson")
            if w.technology in ("onwind", "solar", "solar-hsat")
            else resources("regions_offshore_base_s_{clusters}.geojson")
        ),
        cutout=lambda w: input_cutout(
            w, config_provider("renewable", w.technology, "cutout")(w)
        ),
    output:
        availability_matrix=resources(
            "availability_matrix_MD-UA_{clusters}_{technology}.nc"
        ),
    log:
        logs("determine_availability_matrix_MD_UA_{clusters}_{technology}.log"),
    benchmark:
        benchmarks("determine_availability_matrix_MD_UA_{clusters}_{technology}")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 5000,
    script:
        scripts("determine_availability_matrix_MD_UA.py")


# Optional input when having Ukraine (UA) or Moldova (MD) in the countries list
def input_ua_md_availability_matrix(w):
    countries = set(config_provider("countries")(w))
    if {"UA", "MD"}.intersection(countries):
        return {
            "availability_matrix_MD_UA": resources(
                "availability_matrix_MD-UA_{clusters}_{technology}.nc"
            )
        }
    return {}


rule determine_availability_matrix:
    message:
        "Determining availability matrix for {wildcards.clusters} clusters and {wildcards.technology} technology"
    params:
        renewable=config_provider("renewable"),
    input:
        unpack(input_ua_md_availability_matrix),
        corine=ancient(rules.retrieve_corine.output["tif_file"]),
        natura=lambda w: (
            f"{NATURA_DATASET['folder']}/natura.tiff"
            if config_provider("renewable", w.technology, "natura")(w)
            else []
        ),
        luisa=rules.retrieve_luisa_land_cover.output["tif"],
        gebco=ancient(
            lambda w: (
                rules.retrieve_gebco.output["gebco"]
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
            resources("regions_onshore_base_s_{clusters}.geojson")
            if w.technology in ("onwind", "solar", "solar-hsat")
            else resources("regions_offshore_base_s_{clusters}.geojson")
        ),
        cutout=lambda w: input_cutout(
            w, config_provider("renewable", w.technology, "cutout")(w)
        ),
    output:
        resources("availability_matrix_{clusters}_{technology}.nc"),
    log:
        logs("determine_availability_matrix_{clusters}_{technology}.log"),
    benchmark:
        benchmarks("determine_availability_matrix_{clusters}_{technology}")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 5000,
    script:
        scripts("determine_availability_matrix.py")


rule build_renewable_profiles:
    message:
        "Building renewable profiles for {wildcards.clusters} clusters and {wildcards.technology} technology"
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        renewable=config_provider("renewable"),
    input:
        availability_matrix=resources("availability_matrix_{clusters}_{technology}.nc"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        distance_regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        resource_regions=lambda w: (
            resources("regions_onshore_base_s_{clusters}.geojson")
            if w.technology in ("onwind", "solar", "solar-hsat")
            else resources("regions_offshore_base_s_{clusters}.geojson")
        ),
        cutout=lambda w: input_cutout(
            w, config_provider("renewable", w.technology, "cutout")(w)
        ),
    output:
        profile=resources("profile_{clusters}_{technology}.nc"),
        class_regions=resources("regions_by_class_{clusters}_{technology}.geojson"),
    log:
        logs("build_renewable_profile_{clusters}_{technology}.log"),
    benchmark:
        benchmarks("build_renewable_profile_{clusters}_{technology}")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 5000,
    wildcard_constraints:
        technology="(?!hydro).*",  # Any technology other than hydro
    script:
        scripts("build_renewable_profiles.py")


rule build_monthly_prices:
    message:
        "Building monthly fuel and CO2 prices"
    input:
        co2_price_raw="data/validation/emission-spot-primary-market-auction-report-2019-data.xls",
        fuel_price_raw="data/validation/energy-price-trends-xlsx-5619002.xlsx",
    output:
        co2_price=resources("co2_price.csv"),
        fuel_price=resources("monthly_fuel_price.csv"),
    log:
        logs("build_monthly_prices.log"),
    benchmark:
        benchmarks("build_monthly_prices")
    threads: 1
    resources:
        mem_mb=5000,
    script:
        scripts("build_monthly_prices.py")


if COUNTRY_RUNOFF_DATASET["source"] == "build":

    # This rule uses one or multiple cutouts.
    # To update the output files to include a new year, e.g. 2025 using an existing cutout,
    # either create a new cutout covering the whole timespan or add another cutout that covers the additional year(s).
    # E.g. cutouts=[<cutout for 1940-2024>, <cutout for 2025-2025>]
    rule build_country_runoff:
        input:
            cutouts=["cutouts/europe-1940-2024-era5.nc"],
            country_shapes=resources("country_shapes.geojson"),
        output:
            era5_runoff=COUNTRY_RUNOFF_DATASET["folder"] / "era5-runoff-per-country.csv",
        log:
            logs("build_country_runoff.log"),
        benchmark:
            benchmarks("build_country_runoff")
        script:
            scripts("build_country_runoff.py")


rule build_hydro_profile:
    message:
        "Building hydropower profile"
    params:
        hydro=config_provider("renewable", "hydro"),
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        country_shapes=resources("country_shapes.geojson"),
        eia_hydro_generation="data/eia_hydro_annual_generation.csv",
        eia_hydro_capacity="data/eia_hydro_annual_capacity.csv",
        era5_runoff=f"{COUNTRY_RUNOFF_DATASET['folder']}/era5-runoff-per-country.csv",
        cutout=lambda w: input_cutout(
            w, config_provider("renewable", "hydro", "cutout")(w)
        ),
    output:
        profile=resources("profile_hydro.nc"),
    log:
        logs("build_hydro_profile.log"),
    benchmark:
        benchmarks("build_hydro_profile")
    resources:
        mem_mb=5000,
    script:
        scripts("build_hydro_profile.py")


rule build_line_rating:
    message:
        "Building dynamic line ratings"
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        base_network=resources("networks/base.nc"),
        cutout=lambda w: input_cutout(
            w, config_provider("lines", "dynamic_line_rating", "cutout")(w)
        ),
    output:
        output=resources("dlr.nc"),
    log:
        logs("build_line_rating.log"),
    benchmark:
        benchmarks("build_line_rating")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=config["atlite"].get("nprocesses", 4) * 1000,
    script:
        scripts("build_line_rating.py")


rule build_transmission_projects:
    message:
        "Building transmission projects"
    params:
        transmission_projects=config_provider("transmission_projects"),
        line_factor=config_provider("lines", "length_factor"),
        s_max_pu=config_provider("lines", "s_max_pu"),
    input:
        base_network=resources("networks/base.nc"),
        offshore_shapes=resources("offshore_shapes.geojson"),
        europe_shape=resources("europe_shape.geojson"),
        transmission_projects=lambda w: [
            "data/transmission_projects/" + name
            for name, include in config_provider("transmission_projects", "include")(
                w
            ).items()
            if include
        ],
    output:
        new_lines=resources("transmission_projects/new_lines.csv"),
        new_links=resources("transmission_projects/new_links.csv"),
        adjust_lines=resources("transmission_projects/adjust_lines.csv"),
        adjust_links=resources("transmission_projects/adjust_links.csv"),
        new_buses=resources("transmission_projects/new_buses.csv"),
    log:
        logs("build_transmission_projects.log"),
    benchmark:
        benchmarks("build_transmission_projects")
    resources:
        mem_mb=4000,
    threads: 1
    script:
        scripts("build_transmission_projects.py")


rule add_transmission_projects_and_dlr:
    message:
        "Adding transmission projects and dynamic line ratings"
    params:
        transmission_projects=config_provider("transmission_projects"),
        dlr=config_provider("lines", "dynamic_line_rating"),
        s_max_pu=config_provider("lines", "s_max_pu"),
    input:
        network=resources("networks/base.nc"),
        dlr=lambda w: (
            resources("dlr.nc")
            if config_provider("lines", "dynamic_line_rating", "activate")(w)
            else []
        ),
        transmission_projects=lambda w: (
            [
                resources("transmission_projects/new_buses.csv"),
                resources("transmission_projects/new_lines.csv"),
                resources("transmission_projects/new_links.csv"),
                resources("transmission_projects/adjust_lines.csv"),
                resources("transmission_projects/adjust_links.csv"),
            ]
            if config_provider("transmission_projects", "enable")(w)
            else []
        ),
    output:
        network=resources("networks/base_extended.nc"),
    log:
        logs("add_transmission_projects_and_dlr.log"),
    benchmark:
        benchmarks("add_transmission_projects_and_dlr")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("add_transmission_projects_and_dlr.py")


def input_class_regions(w):
    return {
        f"class_regions_{tech}": resources(
            f"regions_by_class_{{clusters}}_{tech}.geojson"
        )
        for tech in set(config_provider("electricity", "renewable_carriers")(w))
        - {"hydro"}
    }


rule build_electricity_demand_base:
    message:
        "Building electricity demand time series for base network"
    params:
        distribution_key=config_provider("load", "distribution_key"),
    input:
        base_network=resources("networks/base_s.nc"),
        regions=resources("regions_onshore_base_s.geojson"),
        nuts3=resources("nuts3_shapes.geojson"),
        load=resources("electricity_demand.csv"),
    output:
        resources("electricity_demand_base_s.nc"),
    log:
        logs("build_electricity_demand_base_s.log"),
    benchmark:
        benchmarks("build_electricity_demand_base_s")
    resources:
        mem_mb=5000,
    script:
        scripts("build_electricity_demand_base.py")


rule build_hac_features:
    message:
        "Aggregate all rastered cutout data to base regions Voronoi cells."
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        features=config_provider("clustering", "cluster_network", "hac_features"),
    input:
        cutout=lambda w: input_cutout(w),
        regions=resources("regions_onshore_base_s.geojson"),
    output:
        resources("hac_features.nc"),
    log:
        logs("build_hac_features.log"),
    benchmark:
        benchmarks("build_hac_features")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=10000,
    script:
        scripts("build_hac_features.py")


rule process_cost_data:
    params:
        costs=config_provider("costs"),
        max_hours=config_provider("electricity", "max_hours"),
    input:
        network=resources("networks/base_s.nc"),
        costs=rules.retrieve_cost_data.output["costs"],
        custom_costs=config_provider("costs", "custom_cost_fn"),
    output:
        resources("costs_{planning_horizons}_processed.csv"),
    log:
        logs("build_cost_data_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_cost_data_{planning_horizons}")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("process_cost_data.py")


rule simplify_network:
    message:
        "Simplifying network"
    params:
        countries=config_provider("countries"),
        mode=config_provider("clustering", "mode"),
        administrative=config_provider("clustering", "administrative"),
        simplify_network=config_provider("clustering", "simplify_network"),
        cluster_network=config_provider("clustering", "cluster_network"),
        aggregation_strategies=config_provider(
            "clustering", "aggregation_strategies", default={}
        ),
        p_max_pu=config_provider("links", "p_max_pu", default=1.0),
        p_min_pu=config_provider("links", "p_min_pu", default=-1.0),
    input:
        network=resources("networks/base_extended.nc"),
        regions_onshore=resources("regions_onshore.geojson"),
        regions_offshore=resources("regions_offshore.geojson"),
        admin_shapes=resources("admin_shapes.geojson"),
    output:
        network=resources("networks/base_s.nc"),
        regions_onshore=resources("regions_onshore_base_s.geojson"),
        regions_offshore=resources("regions_offshore_base_s.geojson"),
        busmap=resources("busmap_base_s.csv"),
    log:
        logs("simplify_network.log"),
    benchmark:
        benchmarks("simplify_network_b")
    threads: 1
    resources:
        mem_mb=12000,
    script:
        scripts("simplify_network.py")


# Optional input when using custom busmaps - Needs to be tailored to selected base_network
def input_custom_busmap(w):

    custom_busmap = []
    custom_busshapes = []

    mode = config_provider("clustering", "mode", default="busmap")(w)

    if mode == "custom_busmap":
        base_network = config_provider("electricity", "base_network")(w)
        custom_busmap = f"data/busmaps/base_s_{w.clusters}_{base_network}.csv"

    if mode == "custom_busshapes":
        base_network = config_provider("electricity", "base_network")(w)
        custom_busshapes = f"data/busshapes/base_s_{w.clusters}_{base_network}.geojson"

    return {
        "custom_busmap": custom_busmap,
        "custom_busshapes": custom_busshapes,
    }


rule cluster_network:
    message:
        "Clustering network to {wildcards.clusters} clusters"
    params:
        countries=config_provider("countries"),
        mode=config_provider("clustering", "mode"),
        administrative=config_provider("clustering", "administrative"),
        cluster_network=config_provider("clustering", "cluster_network"),
        aggregation_strategies=config_provider(
            "clustering", "aggregation_strategies", default={}
        ),
        focus_weights=config_provider("clustering", "focus_weights", default=None),
        renewable_carriers=config_provider("electricity", "renewable_carriers"),
        conventional_carriers=config_provider(
            "electricity", "conventional_carriers", default=[]
        ),
        max_hours=config_provider("electricity", "max_hours"),
        length_factor=config_provider("lines", "length_factor"),
        cluster_mode=config_provider("clustering", "mode"),
        copperplate_regions=config_provider("clustering", "copperplate_regions"),
    input:
        unpack(input_custom_busmap),
        network=resources("networks/base_s.nc"),
        admin_shapes=resources("admin_shapes.geojson"),
        bidding_zones=lambda w: (
            resources("bidding_zones.geojson")
            if config_provider("clustering", "mode")(w) == "administrative"
            else []
        ),
        regions_onshore=resources("regions_onshore_base_s.geojson"),
        regions_offshore=resources("regions_offshore_base_s.geojson"),
        hac_features=lambda w: (
            resources("hac_features.nc")
            if config_provider("clustering", "cluster_network", "algorithm")(w)
            == "hac"
            else []
        ),
        load=resources("electricity_demand_base_s.nc"),
    output:
        network=resources("networks/base_s_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        busmap=resources("busmap_base_s_{clusters}.csv"),
        linemap=resources("linemap_base_s_{clusters}.csv"),
    log:
        logs("cluster_network_base_s_{clusters}.log"),
    benchmark:
        benchmarks("cluster_network_base_s_{clusters}")
    threads: 1
    resources:
        mem_mb=10000,
    script:
        scripts("cluster_network.py")


def input_profile_tech(w):
    return {
        f"profile_{tech}": resources(
            "profile_{clusters}_" + tech + ".nc"
            if tech != "hydro"
            else f"profile_{tech}.nc"
        )
        for tech in config_provider("electricity", "renewable_carriers")(w)
    }


def input_conventional(w):
    carriers = [
        *config_provider("electricity", "conventional_carriers")(w),
        *config_provider("electricity", "extendable_carriers", "Generator")(w),
    ]
    return {
        f"conventional_{carrier}_{attr}": fn
        for carrier, d in config_provider("conventional", default={})(w).items()
        if carrier in carriers
        for attr, fn in d.items()
        if str(fn).startswith("data/")
    }


rule add_electricity:
    message:
        "Adding electricity to network with {wildcards.clusters} clusters"
    params:
        line_length_factor=config_provider("lines", "length_factor"),
        link_length_factor=config_provider("links", "length_factor"),
        scaling_factor=config_provider("load", "scaling_factor"),
        countries=config_provider("countries"),
        snapshots=config_provider("snapshots"),
        renewable=config_provider("renewable"),
        electricity=config_provider("electricity"),
        conventional=config_provider("conventional"),
        foresight=config_provider("foresight"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        consider_efficiency_classes=config_provider(
            "clustering", "consider_efficiency_classes"
        ),
        aggregation_strategies=config_provider("clustering", "aggregation_strategies"),
        exclude_carriers=config_provider("clustering", "exclude_carriers"),
    input:
        unpack(input_profile_tech),
        unpack(input_class_regions),
        unpack(input_conventional),
        base_network=resources("networks/base_s_{clusters}.nc"),
        costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}_processed.csv"
        ),
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        powerplants=resources("powerplants_s_{clusters}.csv"),
        hydro_capacities=ancient("data/hydro_capacities.csv"),
        unit_commitment="data/unit_commitment.csv",
        fuel_price=lambda w: (
            resources("monthly_fuel_price.csv")
            if config_provider("conventional", "dynamic_fuel_price")(w)
            else []
        ),
        load=resources("electricity_demand_base_s.nc"),
        busmap=resources("busmap_base_s_{clusters}.csv"),
    output:
        resources("networks/base_s_{clusters}_elec.nc"),
    log:
        logs("add_electricity_{clusters}.log"),
    benchmark:
        benchmarks("add_electricity_{clusters}")
    threads: 1
    resources:
        mem_mb=10000,
    script:
        scripts("add_electricity.py")


rule prepare_network:
    message:
        "Preparing network for model with {wildcards.clusters} clusters and options {wildcards.opts}"
    params:
        time_resolution=config_provider("clustering", "temporal", "resolution_elec"),
        links=config_provider("links"),
        lines=config_provider("lines"),
        co2base=config_provider("electricity", "co2base"),
        co2limit_enable=config_provider("electricity", "co2limit_enable", default=False),
        co2limit=config_provider("electricity", "co2limit"),
        gaslimit_enable=config_provider("electricity", "gaslimit_enable", default=False),
        gaslimit=config_provider("electricity", "gaslimit"),
        emission_prices=config_provider("costs", "emission_prices"),
        adjustments=config_provider("adjustments", "electricity"),
        autarky=config_provider("electricity", "autarky", default={}),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        transmission_limit=config_provider("electricity", "transmission_limit"),
    input:
        resources("networks/base_s_{clusters}_elec.nc"),
        costs=lambda w: resources(
            f"costs_{config_provider('costs', 'year')(w)}_processed.csv"
        ),
        co2_price=lambda w: resources("co2_price.csv") if "Ept" in w.opts else [],
    output:
        resources("networks/base_s_{clusters}_elec_{opts}.nc"),
    log:
        logs("prepare_network_base_s_{clusters}_elec_{opts}.log"),
    benchmark:
        benchmarks("prepare_network_base_s_{clusters}_elec_{opts}")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("prepare_network.py")


if (
    config["electricity"]["base_network"] == "osm"
    and config["data"]["osm"]["source"] == "build"
):

    rule clean_osm_data:
        message:
            "Cleaning raw OSM data for countries: " + ", ".join(config["countries"])
        input:
            cables_way=expand(
                f"{OSM_DATASET['folder']}/{{country}}/cables_way.json",
                country=config_provider("countries"),
            ),
            lines_way=expand(
                f"{OSM_DATASET['folder']}/{{country}}/lines_way.json",
                country=config_provider("countries"),
            ),
            routes_relation=expand(
                f"{OSM_DATASET['folder']}/{{country}}/routes_relation.json",
                country=config_provider("countries"),
            ),
            substations_way=expand(
                f"{OSM_DATASET['folder']}/{{country}}/substations_way.json",
                country=config_provider("countries"),
            ),
            substations_relation=expand(
                f"{OSM_DATASET['folder']}/{{country}}/substations_relation.json",
                country=config_provider("countries"),
            ),
            offshore_shapes=resources("offshore_shapes.geojson"),
            country_shapes=resources("country_shapes.geojson"),
        output:
            substations=resources(f"osm/clean/substations.geojson"),
            substations_polygon=resources(f"osm/clean/substations_polygon.geojson"),
            converters_polygon=resources(f"osm/clean/converters_polygon.geojson"),
            lines=resources(f"osm/clean/lines.geojson"),
            links=resources(f"osm/clean/links.geojson"),
        log:
            logs("clean_osm_data.log"),
        benchmark:
            benchmarks("clean_osm_data")
        threads: 1
        resources:
            mem_mb=4000,
        script:
            scripts("clean_osm_data.py")

    rule build_osm_network:
        message:
            "Building OSM network"
        params:
            countries=config_provider("countries"),
            voltages=config_provider("electricity", "voltages"),
            line_types=config_provider("lines", "types"),
        input:
            substations=resources(f"osm/clean/substations.geojson"),
            substations_polygon=resources(f"osm/clean/substations_polygon.geojson"),
            converters_polygon=resources(f"osm/clean/converters_polygon.geojson"),
            lines=resources(f"osm/clean/lines.geojson"),
            links=resources(f"osm/clean/links.geojson"),
            country_shapes=resources("country_shapes.geojson"),
        output:
            lines=resources(f"osm/build/lines.csv"),
            links=resources(f"osm/build/links.csv"),
            converters=resources(f"osm/build/converters.csv"),
            transformers=resources(f"osm/build/transformers.csv"),
            substations=resources(f"osm/build/buses.csv"),
            lines_geojson=resources(f"osm/geojson/lines.geojson"),
            links_geojson=resources(f"osm/geojson/links.geojson"),
            converters_geojson=resources(f"osm/geojson/converters.geojson"),
            transformers_geojson=resources(f"osm/geojson/transformers.geojson"),
            substations_geojson=resources(f"osm/geojson/buses.geojson"),
            stations_polygon=resources(f"osm/geojson/stations_polygon.geojson"),
            buses_polygon=resources(f"osm/geojson/buses_polygon.geojson"),
        log:
            logs("build_osm_network.log"),
        benchmark:
            benchmarks("build_osm_network")
        threads: 1
        resources:
            mem_mb=4000,
        script:
            scripts("build_osm_network.py")


if config["electricity"]["base_network"] == "tyndp":

    rule build_tyndp_network:
        message:
            "Building TYNDP network"
        params:
            countries=config_provider("countries"),
        input:
            reference_grid=rules.retrieve_tyndp.output.reference_grid,
            buses=rules.retrieve_tyndp.output.nodes,
            bidding_shapes=resources("bidding_zones.geojson"),
        output:
            lines=resources("tyndp/build/lines.csv"),
            links=resources("tyndp/build/links.csv"),
            converters=resources("tyndp/build/converters.csv"),
            transformers=resources("tyndp/build/transformers.csv"),
            substations=resources("tyndp/build/buses.csv"),
            substations_h2=resources("tyndp/build/buses_h2.csv"),
            lines_geojson=resources("tyndp/build/geojson/lines.geojson"),
            links_geojson=resources("tyndp/build/geojson/links.geojson"),
            converters_geojson=resources("tyndp/build/geojson/converters.geojson"),
            transformers_geojson=resources("tyndp/build/geojson/transformers.geojson"),
            substations_geojson=resources("tyndp/build/geojson/buses.geojson"),
            substations_h2_geojson=resources("tyndp/build/geojson/buses_h2.geojson"),
        log:
            logs("build_tyndp_network.log"),
        benchmark:
            benchmarks("build_tyndp_network")
        threads: 1
        resources:
            mem_mb=4000,
        script:
            scripts("build_tyndp_network.py")
