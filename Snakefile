
from os.path import exists
from shutil import copyfile

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")

configfile: "config.yaml"


wildcard_constraints:
    lv="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*",
    clusters="[0-9]+m?",
    opts="[-+a-zA-Z0-9]*",
    sector_opts="[-+a-zA-Z0-9\.\s]*"


SDIR = config['summary_dir'] + '/' + config['run']
RDIR = config['results_dir'] + config['run']
CDIR = config['costs_dir']


subworkflow pypsaeur:
    workdir: "../pypsa-eur"
    snakefile: "../pypsa-eur/Snakefile"
    configfile: "../pypsa-eur/config.yaml"

rule all:
    input: SDIR + '/graphs/costs.pdf'


rule solve_all_networks:
    input:
        expand(RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
               **config['scenario'])


rule prepare_sector_networks:
    input:
        expand(RDIR + "/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
               **config['scenario'])

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

if config.get('retrieve_sector_databundle', True):
    rule retrieve_sector_databundle:
        output: *datafiles
        log: "logs/retrieve_sector_databundle.log"
        script: 'scripts/retrieve_sector_databundle.py'


rule build_population_layouts:
    input:
        nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson'),
        urban_percent="data/urban_percent.csv"
    output:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc"
    resources: mem_mb=20000
    benchmark: "benchmarks/build_population_layouts"
    threads: 8
    script: "scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur('resources/regions_onshore_elec_s{simpl}_{clusters}.geojson')
    output:
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv"
    resources: mem_mb=10000
    benchmark: "benchmarks/build_clustered_population_layouts/s{simpl}_{clusters}"
    script: "scripts/build_clustered_population_layouts.py"


rule build_simplified_population_layouts:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur('resources/regions_onshore_elec_s{simpl}.geojson')
    output:
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}.csv"
    resources: mem_mb=10000
    benchmark: "benchmarks/build_clustered_population_layouts/s{simpl}"
    script: "scripts/build_clustered_population_layouts.py"


if config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]:

    datafiles = [
        "IGGIELGN_LNGs.geojson",
        "IGGIELGN_BorderPoints.geojson",
        "IGGIELGN_Productions.geojson",
        "IGGIELGN_PipeSegments.geojson",
    ]


    rule retrieve_gas_infrastructure_data:
        output: expand("data/gas_network/scigrid-gas/data/{files}", files=datafiles)
        script: 'scripts/retrieve_gas_infrastructure_data.py'


    rule build_gas_network:
        input:
            gas_network="data/gas_network/scigrid-gas/data/IGGIELGN_PipeSegments.geojson"
        output:
            cleaned_gas_network="resources/gas_network.csv"
        resources: mem_mb=4000
        script: "scripts/build_gas_network.py"


    rule build_gas_input_locations:
        input:
            lng="data/gas_network/scigrid-gas/data/IGGIELGN_LNGs.geojson",
            entry="data/gas_network/scigrid-gas/data/IGGIELGN_BorderPoints.geojson",
            production="data/gas_network/scigrid-gas/data/IGGIELGN_Productions.geojson",
            planned_lng="data/gas_network/planned_LNGs.csv",
            regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"),
            regions_offshore=pypsaeur('resources/regions_offshore_elec_s{simpl}_{clusters}.geojson')
        output:
            gas_input_nodes="resources/gas_input_locations_s{simpl}_{clusters}.geojson",
            gas_input_nodes_simplified="resources/gas_input_locations_s{simpl}_{clusters}_simplified.csv"
        resources: mem_mb=2000,
        script: "scripts/build_gas_input_locations.py"


    rule cluster_gas_network:
        input:
            cleaned_gas_network="resources/gas_network.csv",
            regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"),
            regions_offshore=pypsaeur("resources/regions_offshore_elec_s{simpl}_{clusters}.geojson")
        output:
            clustered_gas_network="resources/gas_network_elec_s{simpl}_{clusters}.csv"
        resources: mem_mb=4000
        script: "scripts/cluster_gas_network.py"


    gas_infrastructure = {**rules.cluster_gas_network.output, **rules.build_gas_input_locations.output}
else:
    gas_infrastructure = {}


rule build_heat_demands:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson")
    output:
        heat_demand_urban="resources/heat_demand_urban_elec_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/heat_demand_rural_elec_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/heat_demand_total_elec_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    benchmark: "benchmarks/build_heat_demands/s{simpl}_{clusters}"
    script: "scripts/build_heat_demand.py"


rule build_temperature_profiles:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson")
    output:
        temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    benchmark: "benchmarks/build_temperature_profiles/s{simpl}_{clusters}"
    script: "scripts/build_temperature_profiles.py"


rule build_cop_profiles:
    input:
        temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc"
    output:
        cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/cop_air_urban_elec_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    benchmark: "benchmarks/build_cop_profiles/s{simpl}_{clusters}"
    script: "scripts/build_cop_profiles.py"


rule build_solar_thermal_profiles:
    input:
        pop_layout_total="resources/pop_layout_total.nc",
        pop_layout_urban="resources/pop_layout_urban.nc",
        pop_layout_rural="resources/pop_layout_rural.nc",
        regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson")
    output:
        solar_thermal_total="resources/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/solar_thermal_urban_elec_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/solar_thermal_rural_elec_s{simpl}_{clusters}.nc"
    resources: mem_mb=20000
    benchmark: "benchmarks/build_solar_thermal_profiles/s{simpl}_{clusters}"
    script: "scripts/build_solar_thermal_profiles.py"


def input_eurostat(w):
    # 2016 includes BA, 2017 does not
    report_year = config["energy"]["eurostat_report_year"]
    return f"data/eurostat-energy_balances-june_{report_year}_edition"

rule build_energy_totals:
    input:
        nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson'),
        co2="data/eea/UNFCCC_v23.csv",
        swiss="data/switzerland-sfoe/switzerland-new_format.csv",
        idees="data/jrc-idees-2015",
        district_heat_share='data/district_heat_share.csv',
        eurostat=input_eurostat
    output:
        energy_name='resources/energy_totals.csv',
	    co2_name='resources/co2_totals.csv',
	    transport_name='resources/transport_data.csv'
    threads: 16
    resources: mem_mb=10000
    benchmark: "benchmarks/build_energy_totals"
    script: 'scripts/build_energy_totals.py'


rule build_biomass_potentials:
    input:
        enspreso_biomass=HTTP.remote("https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/ENSPRESO/ENSPRESO_BIOMASS.xlsx", keep_local=True),
        nuts2="data/nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson", # https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/#nuts21
        regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        nuts3_population=pypsaeur("data/bundle/nama_10r_3popgdp.tsv.gz"),
        swiss_cantons=pypsaeur("data/bundle/ch_cantons.csv"),
        swiss_population=pypsaeur("data/bundle/je-e-21.03.02.xls"),
        country_shapes=pypsaeur('resources/country_shapes.geojson')
    output:
        biomass_potentials_all='resources/biomass_potentials_all_s{simpl}_{clusters}.csv',
        biomass_potentials='resources/biomass_potentials_s{simpl}_{clusters}.csv'
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_biomass_potentials_s{simpl}_{clusters}"
    script: 'scripts/build_biomass_potentials.py'


if config["sector"]["biomass_transport"]:
    rule build_biomass_transport_costs:
        input:
            transport_cost_data=HTTP.remote("publications.jrc.ec.europa.eu/repository/bitstream/JRC98626/biomass potentials in europe_web rev.pdf", keep_local=True)
        output:
            biomass_transport_costs="resources/biomass_transport_costs.csv",
        threads: 1
        resources: mem_mb=1000
        benchmark: "benchmarks/build_biomass_transport_costs"
        script: 'scripts/build_biomass_transport_costs.py'
    build_biomass_transport_costs_output = rules.build_biomass_transport_costs.output
else:
    build_biomass_transport_costs_output = {}


rule build_salt_cavern_potentials:
    input:
        salt_caverns="data/h2_salt_caverns_GWh_per_sqkm.geojson",
        regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=pypsaeur("resources/regions_offshore_elec_s{simpl}_{clusters}.geojson"),
    output:
        h2_cavern_potential="resources/salt_cavern_potentials_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=2000
    benchmark: "benchmarks/build_salt_cavern_potentials_s{simpl}_{clusters}"
    script: "scripts/build_salt_cavern_potentials.py"


rule build_ammonia_production:
    input:
        usgs="data/myb1-2017-nitro.xls"
    output:
        ammonia_production="resources/ammonia_production.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_ammonia_production"
    script: 'scripts/build_ammonia_production.py'


rule build_industry_sector_ratios:
    input:
        ammonia_production="resources/ammonia_production.csv",
        idees="data/jrc-idees-2015"
    output:
        industry_sector_ratios="resources/industry_sector_ratios.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industry_sector_ratios"
    script: 'scripts/build_industry_sector_ratios.py'


rule build_industrial_production_per_country:
    input:
        ammonia_production="resources/ammonia_production.csv",
        jrc="data/jrc-idees-2015",
        eurostat="data/eurostat-energy_balances-may_2018_edition",
    output:
        industrial_production_per_country="resources/industrial_production_per_country.csv"
    threads: 8
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_production_per_country"
    script: 'scripts/build_industrial_production_per_country.py'


rule build_industrial_production_per_country_tomorrow:
    input:
        industrial_production_per_country="resources/industrial_production_per_country.csv"
    output:
        industrial_production_per_country_tomorrow="resources/industrial_production_per_country_tomorrow_{planning_horizons}.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_production_per_country_tomorrow_{planning_horizons}"
    script: 'scripts/build_industrial_production_per_country_tomorrow.py'


rule build_industrial_distribution_key:
    input:
        regions_onshore=pypsaeur('resources/regions_onshore_elec_s{simpl}_{clusters}.geojson'),
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        hotmaps_industrial_database="data/Industrial_Database.csv",
    output:
        industrial_distribution_key="resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_distribution_key/s{simpl}_{clusters}"
    script: 'scripts/build_industrial_distribution_key.py'


rule build_industrial_production_per_node:
    input:
        industrial_distribution_key="resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        industrial_production_per_country_tomorrow="resources/industrial_production_per_country_tomorrow_{planning_horizons}.csv"
    output:
        industrial_production_per_node="resources/industrial_production_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_production_per_node/s{simpl}_{clusters}_{planning_horizons}"
    script: 'scripts/build_industrial_production_per_node.py'


rule build_industrial_energy_demand_per_node:
    input:
        industry_sector_ratios="resources/industry_sector_ratios.csv",
        industrial_production_per_node="resources/industrial_production_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        industrial_energy_demand_per_node_today="resources/industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv"
    output:
        industrial_energy_demand_per_node="resources/industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_energy_demand_per_node/s{simpl}_{clusters}_{planning_horizons}"
    script: 'scripts/build_industrial_energy_demand_per_node.py'


rule build_industrial_energy_demand_per_country_today:
    input:
        jrc="data/jrc-idees-2015",
        ammonia_production="resources/ammonia_production.csv",
        industrial_production_per_country="resources/industrial_production_per_country.csv"
    output:
        industrial_energy_demand_per_country_today="resources/industrial_energy_demand_per_country_today.csv"
    threads: 8
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_energy_demand_per_country_today"
    script: 'scripts/build_industrial_energy_demand_per_country_today.py'


rule build_industrial_energy_demand_per_node_today:
    input:
        industrial_distribution_key="resources/industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
        industrial_energy_demand_per_country_today="resources/industrial_energy_demand_per_country_today.csv"
    output:
        industrial_energy_demand_per_node_today="resources/industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=1000
    benchmark: "benchmarks/build_industrial_energy_demand_per_node_today/s{simpl}_{clusters}"
    script: 'scripts/build_industrial_energy_demand_per_node_today.py'


if config["sector"]["retrofitting"]["retro_endogen"]:
    rule build_retro_cost:
        input:
            building_stock="data/retro/data_building_stock.csv",
            data_tabula="data/retro/tabula-calculator-calcsetbuilding.csv",
            air_temperature = "resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
            u_values_PL="data/retro/u_values_poland.csv",
            tax_w="data/retro/electricity_taxes_eu.csv",
            construction_index="data/retro/comparative_level_investment.csv",
            floor_area_missing="data/retro/floor_area_missing.csv",
            clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
            cost_germany="data/retro/retro_cost_germany.csv",
            window_assumptions="data/retro/window_assumptions.csv",
        output:
            retro_cost="resources/retro_cost_elec_s{simpl}_{clusters}.csv",
            floor_area="resources/floor_area_elec_s{simpl}_{clusters}.csv"
        resources: mem_mb=1000
        benchmark: "benchmarks/build_retro_cost/s{simpl}_{clusters}"
        script: "scripts/build_retro_cost.py"
    build_retro_cost_output = rules.build_retro_cost.output
else:
    build_retro_cost_output = {}


rule build_population_weighted_energy_totals:
    input:
        energy_totals='resources/energy_totals.csv',
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv"
    output: "resources/pop_weighted_energy_totals_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=2000
    script: "scripts/build_population_weighted_energy_totals.py"


rule build_transport_demand:
    input: 
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        pop_weighted_energy_totals="resources/pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        transport_data='resources/transport_data.csv',
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
    output: 
        transport_demand="resources/transport_demand_s{simpl}_{clusters}.csv",
        transport_data="resources/transport_data_s{simpl}_{clusters}.csv",
        avail_profile="resources/avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile="resources/dsm_profile_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=2000
    script: "scripts/build_transport_demand.py"


rule prepare_sector_network:
    input:
        overrides="data/override_component_attrs",
        network=pypsaeur('networks/elec_s{simpl}_{clusters}_ec_lv{lv}_{opts}.nc'),
        energy_totals_name='resources/energy_totals.csv',
        pop_weighted_energy_totals="resources/pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        transport_demand="resources/transport_demand_s{simpl}_{clusters}.csv",
        transport_data="resources/transport_data_s{simpl}_{clusters}.csv",
        avail_profile="resources/avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile="resources/dsm_profile_s{simpl}_{clusters}.csv",
        co2_totals_name='resources/co2_totals.csv',
        biomass_potentials='resources/biomass_potentials_s{simpl}_{clusters}.csv',
        heat_profile="data/heat_load_profile_BDEW.csv",
        costs=CDIR + "costs_{planning_horizons}.csv",
        profile_offwind_ac=pypsaeur("resources/profile_offwind-ac.nc"),
        profile_offwind_dc=pypsaeur("resources/profile_offwind-dc.nc"),
        h2_cavern="resources/salt_cavern_potentials_s{simpl}_{clusters}.csv",
        busmap_s=pypsaeur("resources/busmap_elec_s{simpl}.csv"),
        busmap=pypsaeur("resources/busmap_elec_s{simpl}_{clusters}.csv"),
        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
        simplified_pop_layout="resources/pop_layout_elec_s{simpl}.csv",
        industrial_demand="resources/industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        heat_demand_urban="resources/heat_demand_urban_elec_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/heat_demand_rural_elec_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/heat_demand_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc",
        cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/cop_air_urban_elec_s{simpl}_{clusters}.nc",
        solar_thermal_total="resources/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/solar_thermal_urban_elec_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/solar_thermal_rural_elec_s{simpl}_{clusters}.nc",
        **build_retro_cost_output,
        **build_biomass_transport_costs_output,
        **gas_infrastructure
    output: RDIR + '/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc'
    threads: 1
    resources: mem_mb=2000
    benchmark: RDIR + "/benchmarks/prepare_network/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}"
    script: "scripts/prepare_sector_network.py"


rule plot_network:
    input:
        overrides="data/override_component_attrs",
        network=RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc"
    output:
        map=RDIR + "/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
        today=RDIR + "/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}-today.pdf"
    threads: 2
    resources: mem_mb=10000
    benchmark: RDIR + "/benchmarks/plot_network/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}"
    script: "scripts/plot_network.py"


rule copy_config:
    output: SDIR + '/configs/config.yaml'
    threads: 1
    resources: mem_mb=1000
    benchmark: SDIR + "/benchmarks/copy_config"
    script: "scripts/copy_config.py"


rule make_summary:
    input:
        overrides="data/override_component_attrs",
        networks=expand(
            RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config['scenario']
        ),
        costs=CDIR + "costs_{}.csv".format(config['scenario']['planning_horizons'][0]),
        plots=expand(
            RDIR + "/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config['scenario']
        )
    output:
        nodal_costs=SDIR + '/csvs/nodal_costs.csv',
        nodal_capacities=SDIR + '/csvs/nodal_capacities.csv',
        nodal_cfs=SDIR + '/csvs/nodal_cfs.csv',
        cfs=SDIR + '/csvs/cfs.csv',
        costs=SDIR + '/csvs/costs.csv',
        capacities=SDIR + '/csvs/capacities.csv',
        curtailment=SDIR + '/csvs/curtailment.csv',
        energy=SDIR + '/csvs/energy.csv',
        supply=SDIR + '/csvs/supply.csv',
        supply_energy=SDIR + '/csvs/supply_energy.csv',
        prices=SDIR + '/csvs/prices.csv',
        weighted_prices=SDIR + '/csvs/weighted_prices.csv',
        market_values=SDIR + '/csvs/market_values.csv',
        price_statistics=SDIR + '/csvs/price_statistics.csv',
        metrics=SDIR + '/csvs/metrics.csv'
    threads: 2
    resources: mem_mb=10000
    benchmark: SDIR + "/benchmarks/make_summary"
    script: "scripts/make_summary.py"


rule plot_summary:
    input:
        costs=SDIR + '/csvs/costs.csv',
        energy=SDIR + '/csvs/energy.csv',
        balances=SDIR + '/csvs/supply_energy.csv'
    output:
        costs=SDIR + '/graphs/costs.pdf',
        energy=SDIR + '/graphs/energy.pdf',
        balances=SDIR + '/graphs/balances-energy.pdf'
    threads: 2
    resources: mem_mb=10000
    benchmark: SDIR + "/benchmarks/plot_summary"
    script: "scripts/plot_summary.py"


if config["foresight"] == "overnight":

    rule solve_network:
        input:
            overrides="data/override_component_attrs",
            network=RDIR + "/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
            costs=CDIR + "costs_{planning_horizons}.csv",
            config=SDIR + '/configs/config.yaml'
        output: RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc"
        shadow: "shallow"
        log:
            solver=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
            python=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_python.log",
            memory=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_memory.log"
        threads: config['solving']['solver'].get('threads', 4)
        resources: mem_mb=config['solving']['mem']
        benchmark: RDIR + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}"
        script: "scripts/solve_network.py"


if config["foresight"] == "myopic":

    rule add_existing_baseyear:
        input:
            overrides="data/override_component_attrs",
            network=RDIR + '/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc',
            powerplants=pypsaeur('resources/powerplants.csv'),
            busmap_s=pypsaeur("resources/busmap_elec_s{simpl}.csv"),
            busmap=pypsaeur("resources/busmap_elec_s{simpl}_{clusters}.csv"),
            clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
            costs=CDIR + "costs_{}.csv".format(config['scenario']['planning_horizons'][0]),
            cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
            cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc",
            existing_heating='data/existing_infrastructure/existing_heating_raw.csv',
            country_codes='data/Country_codes.csv',
            existing_solar='data/existing_infrastructure/solar_capacity_IRENA.csv',
            existing_onwind='data/existing_infrastructure/onwind_capacity_IRENA.csv',
            existing_offwind='data/existing_infrastructure/offwind_capacity_IRENA.csv',
        output: RDIR + '/prenetworks-brownfield/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc'
        wildcard_constraints:
            planning_horizons=config['scenario']['planning_horizons'][0] #only applies to baseyear
        threads: 1
        resources: mem_mb=2000
        benchmark: RDIR + '/benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}'
        script: "scripts/add_existing_baseyear.py"


    def solved_previous_horizon(wildcards):
        planning_horizons = config["scenario"]["planning_horizons"]
        i = planning_horizons.index(int(wildcards.planning_horizons))
        planning_horizon_p = str(planning_horizons[i-1])
        return RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_" + planning_horizon_p + ".nc"


    rule add_brownfield:
        input:
            overrides="data/override_component_attrs",
            network=RDIR + '/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc',
            network_p=solved_previous_horizon, #solved network at previous time step
            costs=CDIR + "costs_{planning_horizons}.csv",
            cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
            cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc"
        output: RDIR + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc"
        threads: 4
        resources: mem_mb=10000
        benchmark: RDIR + '/benchmarks/add_brownfield/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}'
        script: "scripts/add_brownfield.py"


    ruleorder: add_existing_baseyear > add_brownfield


    rule solve_network_myopic:
        input:
            overrides="data/override_component_attrs",
            network=RDIR + "/prenetworks-brownfield/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc",
            costs=CDIR + "costs_{planning_horizons}.csv",
            config=SDIR + '/configs/config.yaml'
        output: RDIR + "/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}.nc"
        shadow: "shallow"
        log:
            solver=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
            python=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_python.log",
            memory=RDIR + "/logs/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}_memory.log"
        threads: 4
        resources: mem_mb=config['solving']['mem']
        benchmark: RDIR + "/benchmarks/solve_network/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_{planning_horizons}"
        script: "scripts/solve_network.py"
