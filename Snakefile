
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


rule test_script:
    input:
        expand("resources/heat_demand_urban_elec_s_{clusters}.nc",
                 **config['scenario'])

rule prepare_sector_networks:
    input:
        expand(config['results_dir'] + config['run'] + "/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc",
                 **config['scenario'])


rule build_population_layouts:
    input:
        nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson'),
        urban_percent="data/urban_percent.csv"
    output:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc"
    script: "scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur('resources/regions_onshore_{network}_s{simpl}_{clusters}.geojson')
    output:
        clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv"
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

rule build_industrial_demand:
    input:
        clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv"
    output:
        industrial_demand="resources/industrial_demand_{network}_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industrial_demand.py'




rule prepare_sector_network:
    input:
        network=pypsaeur('networks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc'),
        energy_totals_name='data/energy_totals.csv',
        co2_totals_name='data/co2_totals.csv',
        transport_name='data/transport_data.csv',
        biomass_potentials='data/biomass_potentials.csv',
        timezone_mappings='data/timezone_mappings.csv',
        heat_profile="data/heat_load_profile_DK_AdamJensen.csv",
        costs="data/costs.csv",
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
    output: config['results_dir']  +  config['run'] + '/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc'
    threads: 1
    resources: mem=1000
    benchmark: "benchmarks/prepare_network/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}"
    script: "scripts/prepare_sector_network.py"
