# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule add_existing_baseyear:
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
        costs=config_provider("costs"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        powerplants=resources("powerplants_s_{clusters}.csv"),
        busmap_s=resources("busmap_base_s.csv"),
        busmap=resources("busmap_base_s_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        costs=lambda w: resources(
            "costs_{}.csv".format(
                config_provider("scenario", "planning_horizons", 0)(w)
            )
        ),
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        existing_heating_distribution=resources(
            "existing_heating_distribution_base_s_{clusters}_{planning_horizons}.csv"
        ),
        heating_efficiencies=resources("heating_efficiencies.csv"),
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
    wildcard_constraints:
        # TODO: The first planning_horizon needs to be aligned across scenarios
        # snakemake does not support passing functions to wildcard_constraints
        # reference: https://github.com/snakemake/snakemake/issues/2703
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs(
            "add_existing_baseyear_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "add_existing_baseyear/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_existing_baseyear.py"


def input_profile_tech_brownfield(w):
    return {
        f"profile_{tech}": resources("profile_{clusters}_" + tech + ".nc")
        for tech in config_provider("electricity", "renewable_carriers")(w)
        if tech != "hydro"
    }


rule add_brownfield:
    params:
        H2_retrofit=config_provider("sector", "H2_retrofit"),
        H2_retrofit_capacity_per_CH4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        threshold_capacity=config_provider("existing_capacities", "threshold_capacity"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        carriers=config_provider("electricity", "renewable_carriers"),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
    input:
        unpack(input_profile_tech_brownfield),
        simplify_busmap=resources("busmap_base_s.csv"),
        cluster_busmap=resources("busmap_base_s_{clusters}.csv"),
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        network_p=solved_previous_horizon,  #solved network at previous time step
        costs=resources("costs_{planning_horizons}.csv"),
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
    threads: 4
    resources:
        mem_mb=10000,
    log:
        logs(
            "add_brownfield_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "add_brownfield/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_brownfield.py"


ruleorder: add_existing_baseyear > add_brownfield


rule solve_sector_network_myopic:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential", default=200
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
        ),
        costs=resources("costs_{planning_horizons}.csv"),
    output:
        network=RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"
