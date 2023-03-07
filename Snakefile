# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from os.path import normpath, exists
from shutil import copyfile, move

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")


configfile: "config.yaml"

COSTS = f"data/costs_{config['costs']['year']}.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""

LOGS = "logs/" + RDIR
BENCHMARKS = "benchmarks/" + RDIR
RESOURCES = "resources/" + RDIR if run.get("shared_resources") else "resources/"
RESULTS = "results/" + RDIR


wildcard_constraints:
    simpl="[a-zA-Z0-9]*",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt)",
    opts="[-+a-zA-Z0-9\.]*",
    sector_opts="[-+a-zA-Z0-9\.\s]*",


rule cluster_all_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}.nc",
            **config["scenario"]
        ),


rule extra_components_all_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc",
            **config["scenario"]
        ),


rule prepare_all_networks:
    input:
        expand(
            RESOURCES
            + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule all:
    input:
        RESULTS + "graphs/costs.pdf",


rule solve_all_elec_networks:
    input:
        expand(
            RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule solve_all_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule plot_all_networks:
    input:
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),


if config["enable"].get("prepare_links_p_nom", False):

    rule prepare_links_p_nom:
        output:
            "data/links_p_nom.csv",
        log:
            LOGS + "prepare_links_p_nom.log",
        threads: 1
        resources:
            mem_mb=1500,
        script:
            "scripts/prepare_links_p_nom.py"


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


if config["enable"].get("retrieve_databundle", True):

    rule retrieve_databundle:
        output:
            expand("data/bundle/{file}", file=datafiles),
        log:
            LOGS + "retrieve_databundle.log",
        resources:
            mem_mb=1000,
        script:
            "scripts/retrieve_databundle.py"


rule retrieve_load_data:
    input:
        HTTP.remote(
            "data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv",
            keep_local=True,
            static=True,
        ),
    output:
        "data/load_raw.csv",
    resources:
        mem_mb=5000,
    run:
        move(input[0], output[0])


rule build_load_data:
    input:
        "data/load_raw.csv",
    output:
        RESOURCES + "load.csv",
    log:
        LOGS + "build_load_data.log",
    resources:
        mem_mb=5000,
    script:
        "scripts/build_load_data.py"


rule build_powerplants:
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
    script:
        "scripts/build_powerplants.py"


rule base_network:
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
    script:
        "scripts/base_network.py"


rule build_shapes:
    input:
        naturalearth="data/bundle/naturalearth/ne_10m_admin_0_countries.shp",
        eez="data/bundle/eez/World_EEZ_v8_2014.shp",
        nuts3="data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp",
        nuts3pop="data/bundle/nama_10r_3popgdp.tsv.gz",
        nuts3gdp="data/bundle/nama_10r_3gdp.tsv.gz",
        ch_cantons="data/bundle/ch_cantons.csv",
        ch_popgdp="data/bundle/je-e-21.03.02.xls",
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
    script:
        "scripts/build_shapes.py"


rule build_bus_regions:
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
    script:
        "scripts/build_bus_regions.py"


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        input:
            regions_onshore=RESOURCES + "regions_onshore.geojson",
            regions_offshore=RESOURCES + "regions_offshore.geojson",
        output:
            "cutouts/" + CDIR + "{cutout}.nc",
        log:
            "logs/" + CDIR + "build_cutout/{cutout}.log",
        benchmark:
            "benchmarks/" + CDIR + "build_cutout_{cutout}"
        threads: ATLITE_NPROCESSES
        resources:
            mem_mb=ATLITE_NPROCESSES * 1000,
        script:
            "scripts/build_cutout.py"


if config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            HTTP.remote(
                "zenodo.org/record/6382570/files/{cutout}.nc",
                keep_local=True,
                static=True,
            ),
        output:
            "cutouts/" + CDIR + "{cutout}.nc",
        log:
            "logs/" + CDIR + "retrieve_cutout_{cutout}.log",
        resources:
            mem_mb=5000,
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
            "logs/retrieve_cost_data_{year}.log",
        resources:
            mem_mb=1000,
        run:
            move(input[0], output[0])


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        input:
            natura="data/bundle/natura/Natura2000_end2015.shp",
            cutouts=expand("cutouts/" + CDIR + "{cutouts}.nc", **config["atlite"]),
        output:
            RESOURCES + "natura.tiff",
        resources:
            mem_mb=5000,
        log:
            LOGS + "build_natura_raster.log",
        script:
            "scripts/build_natura_raster.py"


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
        resources:
            mem_mb=5000,
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
    resources:
        mem_mb=5000,
    run:
        move(input[0], output[0])


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
        RESOURCES + "shipdensity_raster.nc",
    log:
        LOGS + "build_ship_raster.log",
    resources:
        mem_mb=5000,
    benchmark:
        BENCHMARKS + "build_ship_raster"
    script:
        "scripts/build_ship_raster.py"


rule build_renewable_profiles:
    input:
        base_network=RESOURCES + "networks/base.nc",
        corine="data/bundle/corine/g250_clc06_V18_5.tif",
        natura=lambda w: (
            RESOURCES + "natura.tiff"
            if config["renewable"][w.technology]["natura"]
            else []
        ),
        gebco=lambda w: (
            "data/bundle/GEBCO_2014_2D.nc"
            if config["renewable"][w.technology].get("max_depth")
            else []
        ),
        ship_density=lambda w: (
            RESOURCES + "shipdensity_raster.nc"
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
    script:
        "scripts/build_renewable_profiles.py"


rule build_hydro_profile:
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
    script:
        "scripts/build_hydro_profile.py"


rule add_electricity:
    input:
        **{
            f"profile_{tech}": RESOURCES + f"profile_{tech}.nc"
            for tech in config["electricity"]["renewable_carriers"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network=RESOURCES + "networks/base.nc",
        tech_costs=COSTS,
        regions=RESOURCES + "regions_onshore.geojson",
        powerplants=RESOURCES + "powerplants.csv",
        hydro_capacities="data/bundle/hydro_capacities.csv",
        geth_hydro_capacities="data/geth2015_hydro_capacities.csv",
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
        mem_mb=5000,
    script:
        "scripts/add_electricity.py"


rule simplify_network:
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
        mem_mb=4000,
    script:
        "scripts/simplify_network.py"


rule cluster_network:
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
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore=RESOURCES
        + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
        busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
        linemap=RESOURCES + "linemap_elec_s{simpl}_{clusters}.csv",
    log:
        LOGS + "cluster_network/elec_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "cluster_network/elec_s{simpl}_{clusters}"
    threads: 1
    resources:
        mem_mb=6000,
    script:
        "scripts/cluster_network.py"


rule add_extra_components:
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
        mem_mb=3000,
    script:
        "scripts/add_extra_components.py"


rule prepare_network:
    input:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc",
        tech_costs=COSTS,
    output:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        LOGS + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log",
    benchmark:
        (
            BENCHMARKS
            + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "scripts/prepare_network.py"


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule solve_network:
    input:
        RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        solver=normpath(LOGS
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
        ),
        python=LOGS
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        memory=LOGS
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log",
    benchmark:
        BENCHMARKS + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    threads: 4
    resources:
        mem_mb=memory,
    shadow:
        "minimal"
    script:
        "scripts/solve_network.py"


rule solve_operations_network:
    input:
        unprepared=RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc",
        optimized=RESULTS
        + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc",
    log:
        solver=normpath(
LOGS
            + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_solver.log"
        ),
        python=LOGS
        + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_python.log",
        memory=LOGS
        + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_memory.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_operations_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 4
    resources:
        mem_mb=(lambda w: 5000 + 372 * int(w.clusters)),
    shadow:
        "minimal"
    script:
        "scripts/solve_operations_network.py"


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

if config["enable"].get("retrieve_sector_databundle", True):

    rule retrieve_sector_databundle:
        output:
            *datafiles,
        log:
            "logs/retrieve_sector_databundle.log",
        script:
            "scripts/retrieve_sector_databundle.py"


rule build_population_layouts:
    input:
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
        urban_percent="data/urban_percent.csv",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        pop_layout_total=RESOURCES + "pop_layout_total.nc",
        pop_layout_urban=RESOURCES + "pop_layout_urban.nc",
        pop_layout_rural=RESOURCES + "pop_layout_rural.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_population_layouts"
    threads: 8
    script:
        "scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total=RESOURCES + "pop_layout_total.nc",
        pop_layout_urban=RESOURCES + "pop_layout_urban.nc",
        pop_layout_rural=RESOURCES + "pop_layout_rural.nc",
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        clustered_pop_layout=RESOURCES
        + "pop_layout_elec_s{simpl}_{clusters}.csv",
    resources:
        mem_mb=10000,
    benchmark:
        "benchmarks/build_clustered_population_layouts/s{simpl}_{clusters}"
    script:
        "scripts/build_clustered_population_layouts.py"


rule build_simplified_population_layouts:
    input:
        pop_layout_total=RESOURCES + "pop_layout_total.nc",
        pop_layout_urban=RESOURCES + "pop_layout_urban.nc",
        pop_layout_rural=RESOURCES + "pop_layout_rural.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}.csv",
    resources:
        mem_mb=10000,
    benchmark:
        "benchmarks/build_clustered_population_layouts/s{simpl}"
    script:
        "scripts/build_clustered_population_layouts.py"


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
        script:
            "scripts/retrieve_gas_infrastructure_data.py"

    rule build_gas_network:
        input:
            gas_network="data/gas_network/scigrid-gas/data/IGGIELGN_PipeSegments.geojson",
        output:
            cleaned_gas_network=RESOURCES + "gas_network.csv",
        resources:
            mem_mb=4000,
        script:
            "scripts/build_gas_network.py"

    rule build_gas_input_locations:
        input:
            lng=HTTP.remote(
                "https://globalenergymonitor.org/wp-content/uploads/2022/09/Europe-Gas-Tracker-August-2022.xlsx",
                keep_local=True,
            ),
            entry="data/gas_network/scigrid-gas/data/IGGIELGN_BorderPoints.geojson",
            production="data/gas_network/scigrid-gas/data/IGGIELGN_Productions.geojson",
            regions_onshore=RESOURCES
            + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore=RESOURCES
            + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
        output:
            gas_input_nodes=RESOURCES
            + "gas_input_locations_s{simpl}_{clusters}.geojson",
            gas_input_nodes_simplified=RESOURCES
            + "gas_input_locations_s{simpl}_{clusters}_simplified.csv",
        resources:
            mem_mb=2000,
        script:
            "scripts/build_gas_input_locations.py"

    rule cluster_gas_network:
        input:
            cleaned_gas_network=RESOURCES + "gas_network.csv",
            regions_onshore=RESOURCES
            + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore=RESOURCES
            + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
        output:
            clustered_gas_network=RESOURCES
            + "gas_network_elec_s{simpl}_{clusters}.csv",
        resources:
            mem_mb=4000,
        script:
            "scripts/cluster_gas_network.py"

    gas_infrastructure = {
        **rules.cluster_gas_network.output,
        **rules.build_gas_input_locations.output,
    }

if not (config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]):
# this is effecively an `else` statement which is however not liked by snakefmt

    gas_infrastructure = {}


rule build_heat_demands:
    input:
        pop_layout=RESOURCES + "pop_layout_{scope}.nc",
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        heat_demand=RESOURCES
        + "heat_demand_{scope}_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    threads: 8
    benchmark:
        "benchmarks/build_heat_demands/{scope}_s{simpl}_{clusters}"
    script:
        "scripts/build_heat_demand.py"


rule build_temperature_profiles:
    input:
        pop_layout=RESOURCES + "pop_layout_{scope}.nc",
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        temp_soil=RESOURCES + "temp_soil_{scope}_elec_s{simpl}_{clusters}.nc",
        temp_air=RESOURCES + "temp_air_{scope}_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    threads: 8
    benchmark:
        "benchmarks/build_temperature_profiles/{scope}_s{simpl}_{clusters}"
    script:
        "scripts/build_temperature_profiles.py"


rule build_cop_profiles:
    input:
        temp_soil_total=RESOURCES
        + "temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural=RESOURCES
        + "temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban=RESOURCES
        + "temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total=RESOURCES
        + "temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural=RESOURCES
        + "temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban=RESOURCES
        + "temp_air_urban_elec_s{simpl}_{clusters}.nc",
    output:
        cop_soil_total=RESOURCES
        + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural=RESOURCES
        + "cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban=RESOURCES
        + "cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES + "cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural=RESOURCES + "cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban=RESOURCES + "cop_air_urban_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    benchmark:
        "benchmarks/build_cop_profiles/s{simpl}_{clusters}"
    script:
        "scripts/build_cop_profiles.py"


rule build_solar_thermal_profiles:
    input:
        pop_layout=RESOURCES + "pop_layout_{scope}.nc",
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        solar_thermal=RESOURCES
        + "solar_thermal_{scope}_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    threads: 16
    benchmark:
        "benchmarks/build_solar_thermal_profiles/{scope}_s{simpl}_{clusters}"
    script:
        "scripts/build_solar_thermal_profiles.py"


def input_eurostat(w):
    # 2016 includes BA, 2017 does not
    report_year = config["energy"]["eurostat_report_year"]
    return f"data/eurostat-energy_balances-june_{report_year}_edition"


rule build_energy_totals:
    input:
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
        co2="data/eea/UNFCCC_v23.csv",
        swiss="data/switzerland-sfoe/switzerland-new_format.csv",
        idees="data/jrc-idees-2015",
        district_heat_share="data/district_heat_share.csv",
        eurostat=input_eurostat,
    output:
        energy_name=RESOURCES + "energy_totals.csv",
        co2_name=RESOURCES + "co2_totals.csv",
        transport_name=RESOURCES + "transport_data.csv",
    threads: 16
    resources:
        mem_mb=10000,
    benchmark:
        "benchmarks/build_energy_totals"
    script:
        "scripts/build_energy_totals.py"


rule build_biomass_potentials:
    input:
        enspreso_biomass=HTTP.remote(
            "https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/ENSPRESO/ENSPRESO_BIOMASS.xlsx",
            keep_local=True,
        ),
        nuts2="data/nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson",  # https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/#nuts21
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        nuts3_population="data/bundle/nama_10r_3popgdp.tsv.gz",
        swiss_cantons="data/bundle/ch_cantons.csv",
        swiss_population="data/bundle/je-e-21.03.02.xls",
        country_shapes=RESOURCES + "country_shapes.geojson",
    output:
        biomass_potentials_all=RESOURCES
        + "biomass_potentials_all_s{simpl}_{clusters}.csv",
        biomass_potentials=RESOURCES
        + "biomass_potentials_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_biomass_potentials_s{simpl}_{clusters}"
    script:
        "scripts/build_biomass_potentials.py"


if config["sector"]["biomass_transport"]:

    rule build_biomass_transport_costs:
        input:
            transport_cost_data=HTTP.remote(
                "publications.jrc.ec.europa.eu/repository/bitstream/JRC98626/biomass potentials in europe_web rev.pdf",
                keep_local=True,
            ),
        output:
            biomass_transport_costs=RESOURCES + "biomass_transport_costs.csv",
        threads: 1
        resources:
            mem_mb=1000,
        benchmark:
            "benchmarks/build_biomass_transport_costs"
        script:
            "scripts/build_biomass_transport_costs.py"

    build_biomass_transport_costs_output = rules.build_biomass_transport_costs.output

if not config["sector"]["biomass_transport"]:
# this is effecively an `else` statement which is however not liked by snakefmt

    build_biomass_transport_costs_output = {}


if config["sector"]["regional_co2_sequestration_potential"]["enable"]:

    rule build_sequestration_potentials:
        input:
            sequestration_potential=HTTP.remote(
                "https://raw.githubusercontent.com/ericzhou571/Co2Storage/main/resources/complete_map_2020_unit_Mt.geojson",
                keep_local=True,
            ),
            regions_onshore=RESOURCES
            + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore=RESOURCES
            + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
        output:
            sequestration_potential=RESOURCES
            + "co2_sequestration_potential_elec_s{simpl}_{clusters}.csv",
        threads: 1
        resources:
            mem_mb=4000,
        benchmark:
            "benchmarks/build_sequestration_potentials_s{simpl}_{clusters}"
        script:
            "scripts/build_sequestration_potentials.py"

    build_sequestration_potentials_output = rules.build_sequestration_potentials.output


if not config["sector"]["regional_co2_sequestration_potential"]["enable"]:
# this is effecively an `else` statement which is however not liked by snakefmt
    build_sequestration_potentials_output = {}


rule build_salt_cavern_potentials:
    input:
        salt_caverns="data/h2_salt_caverns_GWh_per_sqkm.geojson",
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore=RESOURCES
        + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
    output:
        h2_cavern_potential=RESOURCES
        + "salt_cavern_potentials_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    benchmark:
        "benchmarks/build_salt_cavern_potentials_s{simpl}_{clusters}"
    script:
        "scripts/build_salt_cavern_potentials.py"


rule build_ammonia_production:
    input:
        usgs="data/myb1-2017-nitro.xls",
    output:
        ammonia_production=RESOURCES + "ammonia_production.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_ammonia_production"
    script:
        "scripts/build_ammonia_production.py"


rule build_industry_sector_ratios:
    input:
        ammonia_production=RESOURCES + "ammonia_production.csv",
        idees="data/jrc-idees-2015",
    output:
        industry_sector_ratios=RESOURCES + "industry_sector_ratios.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industry_sector_ratios"
    script:
        "scripts/build_industry_sector_ratios.py"


rule build_industrial_production_per_country:
    input:
        ammonia_production=RESOURCES + "ammonia_production.csv",
        jrc="data/jrc-idees-2015",
        eurostat="data/eurostat-energy_balances-may_2018_edition",
    output:
        industrial_production_per_country=RESOURCES
        + "industrial_production_per_country.csv",
    threads: 8
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_production_per_country"
    script:
        "scripts/build_industrial_production_per_country.py"


rule build_industrial_production_per_country_tomorrow:
    input:
        industrial_production_per_country=RESOURCES
        + "industrial_production_per_country.csv",
    output:
        industrial_production_per_country_tomorrow=RESOURCES
        + "industrial_production_per_country_tomorrow_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_production_per_country_tomorrow_{planning_horizons}"
    script:
        "scripts/build_industrial_production_per_country_tomorrow.py"


rule build_industrial_distribution_key:
    input:
        regions_onshore=RESOURCES
        + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        clustered_pop_layout=RESOURCES
        + "pop_layout_elec_s{simpl}_{clusters}.csv",
        hotmaps_industrial_database="data/Industrial_Database.csv",
    output:
        industrial_distribution_key=RESOURCES
        + "industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_distribution_key/s{simpl}_{clusters}"
    script:
        "scripts/build_industrial_distribution_key.py"


rule build_industrial_production_per_node:
    input:
        industrial_distribution_key=RESOURCES
        + "industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        industrial_production_per_country_tomorrow=RESOURCES
        + "industrial_production_per_country_tomorrow_{planning_horizons}.csv",
    output:
        industrial_production_per_node=RESOURCES
        + "industrial_production_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_production_per_node/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_industrial_production_per_node.py"


rule build_industrial_energy_demand_per_node:
    input:
        industry_sector_ratios=RESOURCES + "industry_sector_ratios.csv",
        industrial_production_per_node=RESOURCES
        + "industrial_production_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        industrial_energy_demand_per_node_today=RESOURCES
        + "industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv",
    output:
        industrial_energy_demand_per_node=RESOURCES
        + "industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_energy_demand_per_node/s{simpl}_{clusters}_{planning_horizons}"
    script:
        "scripts/build_industrial_energy_demand_per_node.py"


rule build_industrial_energy_demand_per_country_today:
    input:
        jrc="data/jrc-idees-2015",
        ammonia_production=RESOURCES + "ammonia_production.csv",
        industrial_production_per_country=RESOURCES
        + "industrial_production_per_country.csv",
    output:
        industrial_energy_demand_per_country_today=RESOURCES
        + "industrial_energy_demand_per_country_today.csv",
    threads: 8
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_energy_demand_per_country_today"
    script:
        "scripts/build_industrial_energy_demand_per_country_today.py"


rule build_industrial_energy_demand_per_node_today:
    input:
        industrial_distribution_key=RESOURCES
        + "industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        industrial_energy_demand_per_country_today=RESOURCES
        + "industrial_energy_demand_per_country_today.csv",
    output:
        industrial_energy_demand_per_node_today=RESOURCES
        + "industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        "benchmarks/build_industrial_energy_demand_per_node_today/s{simpl}_{clusters}"
    script:
        "scripts/build_industrial_energy_demand_per_node_today.py"


if config["sector"]["retrofitting"]["retro_endogen"]:

    rule build_retro_cost:
        input:
            building_stock="data/retro/data_building_stock.csv",
            data_tabula="data/retro/tabula-calculator-calcsetbuilding.csv",
            air_temperature=RESOURCES
            + "temp_air_total_elec_s{simpl}_{clusters}.nc",
            u_values_PL="data/retro/u_values_poland.csv",
            tax_w="data/retro/electricity_taxes_eu.csv",
            construction_index="data/retro/comparative_level_investment.csv",
            floor_area_missing="data/retro/floor_area_missing.csv",
            clustered_pop_layout=RESOURCES
            + "pop_layout_elec_s{simpl}_{clusters}.csv",
            cost_germany="data/retro/retro_cost_germany.csv",
            window_assumptions="data/retro/window_assumptions.csv",
        output:
            retro_cost=RESOURCES + "retro_cost_elec_s{simpl}_{clusters}.csv",
            floor_area=RESOURCES + "floor_area_elec_s{simpl}_{clusters}.csv",
        resources:
            mem_mb=1000,
        benchmark:
            "benchmarks/build_retro_cost/s{simpl}_{clusters}"
        script:
            "scripts/build_retro_cost.py"

    build_retro_cost_output = rules.build_retro_cost.output


if not config["sector"]["retrofitting"]["retro_endogen"]:
# this is effecively an `else` statement which is however not liked by snakefmt
    build_retro_cost_output = {}


rule build_population_weighted_energy_totals:
    input:
        energy_totals=RESOURCES + "energy_totals.csv",
        clustered_pop_layout=RESOURCES
        + "pop_layout_elec_s{simpl}_{clusters}.csv",
    output:
        RESOURCES + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "scripts/build_population_weighted_energy_totals.py"


rule build_shipping_demand:
    input:
        ports="data/attributed_ports.json",
        scope=RESOURCES + "europe_shape.geojson",
        regions=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        demand=RESOURCES + "energy_totals.csv",
    output:
        RESOURCES + "shipping_demand_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "scripts/build_shipping_demand.py"


rule build_transport_demand:
    input:
        clustered_pop_layout=RESOURCES
        + "pop_layout_elec_s{simpl}_{clusters}.csv",
        pop_weighted_energy_totals=RESOURCES
        + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        transport_data=RESOURCES + "transport_data.csv",
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        temp_air_total=RESOURCES
        + "temp_air_total_elec_s{simpl}_{clusters}.nc",
    output:
        transport_demand=RESOURCES
        + "transport_demand_s{simpl}_{clusters}.csv",
        transport_data=RESOURCES + "transport_data_s{simpl}_{clusters}.csv",
        avail_profile=RESOURCES + "avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile=RESOURCES + "dsm_profile_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "scripts/build_transport_demand.py"


rule prepare_sector_network:
    params:
        RDIR=RDIR,
    input:
        **build_retro_cost_output,
        **build_biomass_transport_costs_output,
        **gas_infrastructure,
        **build_sequestration_potentials_output,
        overrides="data/override_component_attrs",
        network=RESOURCES
        + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        energy_totals_name=RESOURCES + "energy_totals.csv",
        eurostat=input_eurostat,
        pop_weighted_energy_totals=RESOURCES
        + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        shipping_demand=RESOURCES + "shipping_demand_s{simpl}_{clusters}.csv",
        transport_demand=RESOURCES
        + "transport_demand_s{simpl}_{clusters}.csv",
        transport_data=RESOURCES + "transport_data_s{simpl}_{clusters}.csv",
        avail_profile=RESOURCES + "avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile=RESOURCES + "dsm_profile_s{simpl}_{clusters}.csv",
        co2_totals_name=RESOURCES + "co2_totals.csv",
        co2="data/eea/UNFCCC_v23.csv",
        biomass_potentials=RESOURCES
        + "biomass_potentials_s{simpl}_{clusters}.csv",
        heat_profile="data/heat_load_profile_BDEW.csv",
        costs="data/costs_{}.csv".format(config["costs"]["year"])
        if config["foresight"] == "overnight"
        else "data/costs_{planning_horizons}.csv",
        profile_offwind_ac=RESOURCES + "profile_offwind-ac.nc",
        profile_offwind_dc=RESOURCES + "profile_offwind-dc.nc",
        h2_cavern=RESOURCES + "salt_cavern_potentials_s{simpl}_{clusters}.csv",
        busmap_s=RESOURCES + "busmap_elec_s{simpl}.csv",
        busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
        clustered_pop_layout=RESOURCES
        + "pop_layout_elec_s{simpl}_{clusters}.csv",
        simplified_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}.csv",
        industrial_demand=RESOURCES
        + "industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        heat_demand_urban=RESOURCES
        + "heat_demand_urban_elec_s{simpl}_{clusters}.nc",
        heat_demand_rural=RESOURCES
        + "heat_demand_rural_elec_s{simpl}_{clusters}.nc",
        heat_demand_total=RESOURCES
        + "heat_demand_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_total=RESOURCES
        + "temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural=RESOURCES
        + "temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban=RESOURCES
        + "temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total=RESOURCES
        + "temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural=RESOURCES
        + "temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban=RESOURCES
        + "temp_air_urban_elec_s{simpl}_{clusters}.nc",
        cop_soil_total=RESOURCES
        + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural=RESOURCES
        + "cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban=RESOURCES
        + "cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES + "cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural=RESOURCES + "cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban=RESOURCES + "cop_air_urban_elec_s{simpl}_{clusters}.nc",
        solar_thermal_total=RESOURCES
        + "solar_thermal_total_elec_s{simpl}_{clusters}.nc"
        if config["sector"]["solar_thermal"]
        else [],
        solar_thermal_urban=RESOURCES
        + "solar_thermal_urban_elec_s{simpl}_{clusters}.nc"
        if config["sector"]["solar_thermal"]
        else [],
        solar_thermal_rural=RESOURCES
        + "solar_thermal_rural_elec_s{simpl}_{clusters}.nc"
        if config["sector"]["solar_thermal"]
        else [],
    output:
        RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    benchmark:
        (
            RESULTS
            + "benchmarks/prepare_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "scripts/prepare_sector_network.py"


rule plot_network:
    input:
        overrides="data/override_component_attrs",
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
    output:
        map=RESULTS
        + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
        today=RESULTS
        + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}-today.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RESULTS
            + "benchmarks/plot_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "scripts/plot_network.py"


rule copy_config:
    params:
        RDIR=RDIR,
    output:
        RESULTS + "configs/config.yaml",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        RESULTS + "benchmarks/copy_config"
    script:
        "scripts/copy_config.py"


rule copy_conda_env:
    output:
        RESULTS + "configs/environment.yaml",
    threads: 1
    resources:
        mem_mb=500,
    benchmark:
        RESULTS + "benchmarks/copy_conda_env"
    shell:
        "conda env export -f {output} --no-builds"


rule make_summary:
    params:
        RDIR=RDIR,
    input:
        overrides="data/override_component_attrs",
        networks=expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),
        costs="data/costs_{}.csv".format(config["costs"]["year"])
        if config["foresight"] == "overnight"
        else "data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        plots=expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),
    output:
        nodal_costs=RESULTS + "csvs/nodal_costs.csv",
        nodal_capacities=RESULTS + "csvs/nodal_capacities.csv",
        nodal_cfs=RESULTS + "csvs/nodal_cfs.csv",
        cfs=RESULTS + "csvs/cfs.csv",
        costs=RESULTS + "csvs/costs.csv",
        capacities=RESULTS + "csvs/capacities.csv",
        curtailment=RESULTS + "csvs/curtailment.csv",
        energy=RESULTS + "csvs/energy.csv",
        supply=RESULTS + "csvs/supply.csv",
        supply_energy=RESULTS + "csvs/supply_energy.csv",
        prices=RESULTS + "csvs/prices.csv",
        weighted_prices=RESULTS + "csvs/weighted_prices.csv",
        market_values=RESULTS + "csvs/market_values.csv",
        price_statistics=RESULTS + "csvs/price_statistics.csv",
        metrics=RESULTS + "csvs/metrics.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        RESULTS + "benchmarks/make_summary"
    script:
        "scripts/make_summary.py"


rule plot_summary:
    params:
        RDIR=RDIR,
    input:
        costs=RESULTS + "csvs/costs.csv",
        energy=RESULTS + "csvs/energy.csv",
        balances=RESULTS + "csvs/supply_energy.csv",
        eurostat=input_eurostat,
        country_codes="data/Country_codes.csv",
    output:
        costs=RESULTS + "graphs/costs.pdf",
        energy=RESULTS + "graphs/energy.pdf",
        balances=RESULTS + "graphs/balances-energy.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        RESULTS + "benchmarks/plot_summary"
    script:
        "scripts/plot_summary.py"


if config["foresight"] == "overnight":

    rule solve_sector_network:
        input:
            overrides="data/override_component_attrs",
            network=RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            costs="data/costs_{}.csv".format(config["costs"]["year"]),
            config=RESULTS + "configs/config.yaml",
            #env=RDIR + 'configs/environment.yaml',
        output:
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        shadow:
            "shallow"
        log:
            solver=LOGS
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
            python=LOGS
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
            memory=LOGS
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        threads: config["solving"]["solver"].get("threads", 4)
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RESULTS
                + "benchmarks/solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        script:
            "scripts/solve_sector_network.py"


if config["foresight"] == "myopic":

    rule add_existing_baseyear:
        input:
            overrides="data/override_component_attrs",
            network=RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            powerplants=RESOURCES + "powerplants.csv",
            busmap_s=RESOURCES + "busmap_elec_s{simpl}.csv",
            busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
            clustered_pop_layout=RESOURCES
            + "pop_layout_elec_s{simpl}_{clusters}.csv",
            costs="data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
            cop_soil_total=RESOURCES
            + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
            cop_air_total=RESOURCES
            + "cop_air_total_elec_s{simpl}_{clusters}.nc",
            existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
            country_codes="data/Country_codes.csv",
            existing_solar="data/existing_infrastructure/solar_capacity_IRENA.csv",
            existing_onwind="data/existing_infrastructure/onwind_capacity_IRENA.csv",
            existing_offwind="data/existing_infrastructure/offwind_capacity_IRENA.csv",
        output:
            RESULTS
            + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        wildcard_constraints:
            planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
        threads: 1
        resources:
            mem_mb=2000,
        benchmark:
            (
                RESULTS
                + "benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        script:
            "scripts/add_existing_baseyear.py"

    def solved_previous_horizon(wildcards):
        planning_horizons = config["scenario"]["planning_horizons"]
        i = planning_horizons.index(int(wildcards.planning_horizons))
        planning_horizon_p = str(planning_horizons[i - 1])
        return (
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_"
            + planning_horizon_p
            + ".nc"
        )

    rule add_brownfield:
        input:
            overrides="data/override_component_attrs",
            network=RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            network_p=solved_previous_horizon,  #solved network at previous time step
            costs="data/costs_{planning_horizons}.csv",
            cop_soil_total=RESOURCES
            + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
            cop_air_total=RESOURCES
            + "cop_air_total_elec_s{simpl}_{clusters}.nc",
        output:
            RESULTS
            + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        threads: 4
        resources:
            mem_mb=10000,
        benchmark:
            (
                RESULTS
                + "benchmarks/add_brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        script:
            "scripts/add_brownfield.py"

    ruleorder: add_existing_baseyear > add_brownfield

    rule solve_sector_network_myopic:
        input:
            overrides="data/override_component_attrs",
            network=RESULTS
            + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            costs="data/costs_{planning_horizons}.csv",
            config=RESULTS + "configs/config.yaml",
        output:
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        shadow:
            "shallow"
        log:
            solver=LOGS
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
            python=LOGS
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
            memory=LOGS
            + "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        threads: 4
        resources:
            mem_mb=config["solving"]["mem"],
        benchmark:
            (
                RESULTS
                + "benchmarks/solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
            )
        script:
            "scripts/solve_sector_network.py"
