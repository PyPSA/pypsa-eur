
configfile: "config.yaml"

wildcard_constraints:
    lv="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*",
    clusters="[0-9]+m?",
    sectors="[+a-zA-Z0-9]+",
    opts="[-+a-zA-Z0-9]*",
    sector_opts="[-+a-zA-Z0-9]*"



subworkflow pypsaeur:
    workdir: "../pypsa-eur"
    snakefile: "../pypsa-eur/Snakefile"
    configfile: "../pypsa-eur/config.yaml"

rule all:
    input:
       config['summary_dir'] + '/' + config['run'] + '/graphs/costs.pdf'



rule solve_all_networks:
    input:
        expand(config['results_dir'] + config['run'] + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc",
               **config['scenario'])

rule test_script:
    input:
        expand("resources/heat_demand_urban_elec_s_{clusters}.nc",
                 **config['scenario'])

rule prepare_sector_networks:
    input:
        expand(config['results_dir'] + config['run'] + "/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc",
                 **config['scenario'])


rule build_population_layouts:
    input:
        nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson'),
        urban_percent="data/urban_percent.csv"
    output:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc"
    resources: mem_mb=20000
    script: "scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur('resources/regions_onshore_{network}_s{simpl}_{clusters}.geojson')
    output:
        clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv"
    resources: mem_mb=10000
    script: "scripts/build_clustered_population_layouts.py"


rule build_heat_demands:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur("resources/regions_onshore_{network}_s{simpl}_{clusters}.geojson")
    output:
        heat_demand_urban="resources/heat_demand_urban_{network}_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/heat_demand_rural_{network}_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/heat_demand_total_{network}_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    script: "scripts/build_heat_demand.py"

rule build_temperature_profiles:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur("resources/regions_onshore_{network}_s{simpl}_{clusters}.geojson")
    output:
        temp_soil_total="resources/temp_soil_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_{network}_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_{network}_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_{network}_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_{network}_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_{network}_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    script: "scripts/build_temperature_profiles.py"


rule build_cop_profiles:
    input:
        temp_soil_total="resources/temp_soil_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_{network}_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_{network}_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_{network}_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_{network}_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_{network}_s{simpl}_{clusters}.nc"
    output:
        cop_soil_total="resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/cop_soil_rural_{network}_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/cop_soil_urban_{network}_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cop_air_total_{network}_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/cop_air_rural_{network}_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/cop_air_urban_{network}_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    script: "scripts/build_cop_profiles.py"


rule build_solar_thermal_profiles:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur("resources/regions_onshore_{network}_s{simpl}_{clusters}.geojson")
    output:
        solar_thermal_total="resources/solar_thermal_total_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/solar_thermal_urban_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/solar_thermal_rural_{network}_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    script: "scripts/build_solar_thermal_profiles.py"



rule build_energy_totals:
    input:
        nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson')
    output:
        energy_name='data/energy_totals.csv',
	co2_name='data/co2_totals.csv',
	transport_name='data/transport_data.csv'
    threads: 1
    resources: mem_mb=10000
    script: 'scripts/build_energy_totals.py'

rule build_biomass_potentials:
    input:
        jrc_potentials="data/biomass/JRC Biomass Potentials.xlsx"
    output:
        biomass_potentials='data/biomass_potentials.csv'
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_biomass_potentials.py'


rule build_industry_sector_ratios:
    output:
        industry_sector_ratios="resources/industry_sector_ratios.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industry_sector_ratios.py'


rule build_industrial_demand_per_country:
    input:
        industry_sector_ratios="resources/industry_sector_ratios.csv"
    output:
        industrial_demand_per_country="resources/industrial_demand_per_country.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industrial_demand_per_country.py'


rule build_industrial_demand:
    input:
        clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv",
        industrial_demand_per_country="resources/industrial_demand_per_country.csv"
    output:
        industrial_demand="resources/industrial_demand_{network}_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industrial_demand.py'


rule prepare_sector_network:
    input:
        network=pypsaeur('networks/{network}_s{simpl}_{clusters}_ec_lv{lv}_{opts}.nc'),
        energy_totals_name='data/energy_totals.csv',
        co2_totals_name='data/co2_totals.csv',
        transport_name='data/transport_data.csv',
        biomass_potentials='data/biomass_potentials.csv',
        timezone_mappings='data/timezone_mappings.csv',
        heat_profile="data/heat_load_profile_BDEW.csv",
        costs=config['costs_dir'] + "costs_{planning_horizons}.csv",
        co2_budget="data/co2_budget.csv",
        profile_offwind_ac=pypsaeur("resources/profile_offwind-ac.nc"),
        profile_offwind_dc=pypsaeur("resources/profile_offwind-dc.nc"),
        clustermaps=pypsaeur('resources/clustermaps_{network}_s{simpl}_{clusters}.h5'),
        clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv",
        industrial_demand="resources/industrial_demand_{network}_s{simpl}_{clusters}.csv",
        heat_demand_urban="resources/heat_demand_urban_{network}_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/heat_demand_rural_{network}_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/heat_demand_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_total="resources/temp_soil_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_{network}_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_{network}_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_{network}_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_{network}_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_{network}_s{simpl}_{clusters}.nc",
        cop_soil_total="resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/cop_soil_rural_{network}_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/cop_soil_urban_{network}_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cop_air_total_{network}_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/cop_air_rural_{network}_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/cop_air_urban_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_total="resources/solar_thermal_total_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/solar_thermal_urban_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/solar_thermal_rural_{network}_s{simpl}_{clusters}.nc"
    output: config['results_dir']  +  config['run'] + '/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc'
    threads: 1
    resources: mem_mb=2000
    benchmark: config['results_dir'] + config['run'] + "/benchmarks/prepare_network/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}"
    script: "scripts/prepare_sector_network.py"



rule plot_network:
    input:
        network=config['results_dir'] + config['run'] + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc"
    output:
        map=config['results_dir'] + config['run'] + "/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-costs-all_{co2_budget_name}_{planning_horizons}.pdf",
        today=config['results_dir'] + config['run'] + "/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}-today.pdf"
    threads: 2
    resources: mem_mb=10000
    script: "scripts/plot_network.py"


rule copy_config:
    output:
        config=config['summary_dir'] + '/' + config['run'] + '/configs/config.yaml'
    threads: 1
    resources: mem_mb=1000
    script:
        'scripts/copy_config.py'


rule make_summary:
    input:
        networks=expand(config['results_dir'] + config['run'] + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc",
                 **config['scenario']),
        costs=config['costs_dir'] + "costs_{}.csv".format(config['scenario']['planning_horizons'][0]),
        plots=expand(config['results_dir'] + config['run'] + "/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-costs-all_{co2_budget_name}_{planning_horizons}.pdf",
              **config['scenario'])
        #heat_demand_name='data/heating/daily_heat_demand.h5'
    output:
        nodal_costs=config['summary_dir'] + '/' + config['run'] + '/csvs/nodal_costs.csv',
        nodal_capacities=config['summary_dir'] + '/' + config['run'] + '/csvs/nodal_capacities.csv',
        nodal_cfs=config['summary_dir'] + '/' + config['run'] + '/csvs/nodal_cfs.csv',
        cfs=config['summary_dir'] + '/' + config['run'] + '/csvs/cfs.csv',
        costs=config['summary_dir'] + '/' + config['run'] + '/csvs/costs.csv',
        capacities=config['summary_dir'] + '/' + config['run'] + '/csvs/capacities.csv',
        curtailment=config['summary_dir'] + '/' + config['run'] + '/csvs/curtailment.csv',
        energy=config['summary_dir'] + '/' + config['run'] + '/csvs/energy.csv',
        supply=config['summary_dir'] + '/' + config['run'] + '/csvs/supply.csv',
        supply_energy=config['summary_dir'] + '/' + config['run'] + '/csvs/supply_energy.csv',
        prices=config['summary_dir'] + '/' + config['run'] + '/csvs/prices.csv',
        weighted_prices=config['summary_dir'] + '/' + config['run'] + '/csvs/weighted_prices.csv',
        market_values=config['summary_dir'] + '/' + config['run'] + '/csvs/market_values.csv',
        price_statistics=config['summary_dir'] + '/' + config['run'] + '/csvs/price_statistics.csv',
        metrics=config['summary_dir'] + '/' + config['run'] + '/csvs/metrics.csv'
    threads: 2
    resources: mem_mb=10000
    script:
        'scripts/make_summary.py'


rule plot_summary:
    input:
        costs=config['summary_dir'] + '/' + config['run'] + '/csvs/costs.csv',
        energy=config['summary_dir'] + '/' + config['run'] + '/csvs/energy.csv',
        balances=config['summary_dir'] + '/' + config['run'] + '/csvs/supply_energy.csv'
    output:
        costs=config['summary_dir'] + '/' + config['run'] + '/graphs/costs.pdf',
        energy=config['summary_dir'] + '/' + config['run'] + '/graphs/energy.pdf',
        balances=config['summary_dir'] + '/' + config['run'] + '/graphs/balances-energy.pdf'
    threads: 2
    resources: mem_mb=10000
    script:
        'scripts/plot_summary.py'

if config["foresight"] == "overnight":

    rule solve_network:
        input:
            network=config['results_dir'] + config['run'] + "/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc",
            costs=config['costs_dir'] + "costs_{planning_horizons}.csv",
            config=config['summary_dir'] + '/' + config['run'] + '/configs/config.yaml'
        output: config['results_dir'] + config['run'] + "/postnetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc"
        shadow: "shallow"
        log:
            solver=config['results_dir'] + config['run'] + "/logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}_solver.log",
            python=config['results_dir'] + config['run'] + "/logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}_python.log",
            memory=config['results_dir'] + config['run'] + "/logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}_memory.log"
        benchmark: config['results_dir'] + config['run'] + "/benchmarks/solve_network/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}"
        threads: 4
        resources: mem_mb=config['solving']['mem']
        # group: "solve" # with group, threads is ignored https://bitbucket.org/snakemake/snakemake/issues/971/group-job-description-does-not-contain
        script: "scripts/solve_network.py"


if config["foresight"] == "myopic":

    rule add_existing_baseyear:
        input:
            network=config['results_dir']  +  config['run'] + '/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc',
            powerplants=pypsaeur('resources/powerplants.csv'),
            clustermaps=pypsaeur('resources/clustermaps_{network}_s{simpl}_{clusters}.h5'),
            clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv",
            costs=config['costs_dir'] + "costs_{}.csv".format(config['scenario']['planning_horizons'][0]),
            cop_soil_total="resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc",
            cop_air_total="resources/cop_air_total_{network}_s{simpl}_{clusters}.nc"
        output: config['results_dir']  +  config['run'] + '/prenetworks-brownfield/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc'
        wildcard_constraints:
            planning_horizons=config['scenario']['planning_horizons'][0] #only applies to baseyear
        threads: 1
        resources: mem_mb=2000
        script: "scripts/add_existing_baseyear.py"

    def process_input(wildcards):
        i = config["scenario"]["planning_horizons"].index(int(wildcards.planning_horizons))
        return config['results_dir'] + config['run'] + "/postnetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_" + str(config["scenario"]["planning_horizons"][i-1]) + ".nc"


    rule add_brownfield:
        input:
            network=config['results_dir']  +  config['run'] + '/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc',
            network_p=process_input, #solved network at previous time step
            costs=config['costs_dir'] + "costs_{planning_horizons}.csv",
            cop_soil_total="resources/cop_soil_total_{network}_s{simpl}_{clusters}.nc",
            cop_air_total="resources/cop_air_total_{network}_s{simpl}_{clusters}.nc"

        output: config['results_dir'] + config['run'] + "/prenetworks-brownfield/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc"
        threads: 4
        resources: mem_mb=2000
        script: "scripts/add_brownfield.py"

    ruleorder: add_existing_baseyear > add_brownfield

    rule solve_network_myopic:
        input:
            network=config['results_dir'] + config['run'] + "/prenetworks-brownfield/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc",
            costs=config['costs_dir'] + "costs_{planning_horizons}.csv",
            config=config['summary_dir'] + '/' + config['run'] + '/configs/config.yaml'
        output: config['results_dir'] + config['run'] + "/postnetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}.nc"
        shadow: "shallow"
        log:
            solver=config['results_dir'] + config['run'] + "/logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}_solver.log",
            python=config['results_dir'] + config['run'] + "/logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}_python.log",
            memory=config['results_dir'] + config['run'] + "/logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}_memory.log"
        benchmark: config['results_dir'] + config['run'] + "/benchmarks/solve_network/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{co2_budget_name}_{planning_horizons}"
        threads: 4
        resources: mem_mb=config['solving']['mem']
        script: "scripts/solve_network.py"
