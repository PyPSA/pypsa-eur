# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule build_population_layouts:
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        urban_percent="data/urban_percent.csv",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
    log:
        logs("build_population_layouts.log"),
    resources:
        mem_mb=20000,
    benchmark:
        benchmarks("build_population_layouts")
    threads: 8
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
    log:
        logs("build_clustered_population_layouts_{simpl}_{clusters}.log"),
    resources:
        mem_mb=10000,
    benchmark:
        benchmarks("build_clustered_population_layouts/s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_clustered_population_layouts.py"


rule build_simplified_population_layouts:
    input:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}.geojson"),
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}.csv"),
    resources:
        mem_mb=10000,
    log:
        logs("build_simplified_population_layouts_{simpl}"),
    benchmark:
        benchmarks("build_simplified_population_layouts/s{simpl}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_clustered_population_layouts.py"


if config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]:

    rule build_gas_network:
        input:
            gas_network="data/gas_network/scigrid-gas/data/IGGIELGN_PipeSegments.geojson",
        output:
            cleaned_gas_network=resources("gas_network.csv"),
        resources:
            mem_mb=4000,
        log:
            logs("build_gas_network.log"),
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_gas_network.py"

    rule build_gas_input_locations:
        input:
            lng=HTTP.remote(
                "https://globalenergymonitor.org/wp-content/uploads/2023/07/Europe-Gas-Tracker-2023-03-v3.xlsx",
                keep_local=True,
            ),
            entry="data/gas_network/scigrid-gas/data/IGGIELGN_BorderPoints.geojson",
            production="data/gas_network/scigrid-gas/data/IGGIELGN_Productions.geojson",
            regions_onshore=resources(
                "regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
            regions_offshore=resources(
                "regions_offshore_elec_s{simpl}_{clusters}.geojson"
            ),
        output:
            gas_input_nodes=resources("gas_input_locations_s{simpl}_{clusters}.geojson"),
            gas_input_nodes_simplified=resources(
                "gas_input_locations_s{simpl}_{clusters}_simplified.csv"
            ),
        resources:
            mem_mb=2000,
        log:
            logs("build_gas_input_locations_s{simpl}_{clusters}.log"),
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_gas_input_locations.py"

    rule cluster_gas_network:
        input:
            cleaned_gas_network=resources("gas_network.csv"),
            regions_onshore=resources(
                "regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
            regions_offshore=resources(
                "regions_offshore_elec_s{simpl}_{clusters}.geojson"
            ),
        output:
            clustered_gas_network=resources("gas_network_elec_s{simpl}_{clusters}.csv"),
        resources:
            mem_mb=4000,
        log:
            logs("cluster_gas_network_s{simpl}_{clusters}.log"),
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/cluster_gas_network.py"

    gas_infrastructure = {
        **rules.cluster_gas_network.output,
        **rules.build_gas_input_locations.output,
    }


if not (config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]):
    # this is effecively an `else` statement which is however not liked by snakefmt

    gas_infrastructure = {}


rule build_heat_demands:
    params:
        snapshots=config_provider("snapshots"),
    input:
        pop_layout=resources("pop_layout_{scope}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        heat_demand=resources("heat_demand_{scope}_elec_s{simpl}_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 8
    log:
        logs("build_heat_demands_{scope}_{simpl}_{clusters}.loc"),
    benchmark:
        benchmarks("build_heat_demands/{scope}_s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_heat_demand.py"


rule build_temperature_profiles:
    params:
        snapshots=config_provider("snapshots"),
    input:
        pop_layout=resources("pop_layout_{scope}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        temp_soil=resources("temp_soil_{scope}_elec_s{simpl}_{clusters}.nc"),
        temp_air=resources("temp_air_{scope}_elec_s{simpl}_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 8
    log:
        logs("build_temperature_profiles_{scope}_{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_temperature_profiles/{scope}_s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_temperature_profiles.py"


rule build_cop_profiles:
    params:
        heat_pump_sink_T=config_provider("sector", "heat_pump_sink_T"),
    input:
        temp_soil_total=resources("temp_soil_total_elec_s{simpl}_{clusters}.nc"),
        temp_soil_rural=resources("temp_soil_rural_elec_s{simpl}_{clusters}.nc"),
        temp_soil_urban=resources("temp_soil_urban_elec_s{simpl}_{clusters}.nc"),
        temp_air_total=resources("temp_air_total_elec_s{simpl}_{clusters}.nc"),
        temp_air_rural=resources("temp_air_rural_elec_s{simpl}_{clusters}.nc"),
        temp_air_urban=resources("temp_air_urban_elec_s{simpl}_{clusters}.nc"),
    output:
        cop_soil_total=resources("cop_soil_total_elec_s{simpl}_{clusters}.nc"),
        cop_soil_rural=resources("cop_soil_rural_elec_s{simpl}_{clusters}.nc"),
        cop_soil_urban=resources("cop_soil_urban_elec_s{simpl}_{clusters}.nc"),
        cop_air_total=resources("cop_air_total_elec_s{simpl}_{clusters}.nc"),
        cop_air_rural=resources("cop_air_rural_elec_s{simpl}_{clusters}.nc"),
        cop_air_urban=resources("cop_air_urban_elec_s{simpl}_{clusters}.nc"),
    resources:
        mem_mb=20000,
    log:
        logs("build_cop_profiles_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_cop_profiles/s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_cop_profiles.py"


rule build_solar_thermal_profiles:
    params:
        snapshots=config_provider("snapshots"),
        solar_thermal=config_provider("solar_thermal"),
    input:
        pop_layout=resources("pop_layout_{scope}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        solar_thermal=resources("solar_thermal_{scope}_elec_s{simpl}_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 16
    log:
        logs("build_solar_thermal_profiles_{scope}_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_solar_thermal_profiles/{scope}_s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_solar_thermal_profiles.py"


rule build_energy_totals:
    params:
        countries=config_provider("countries"),
        energy=config_provider("energy"),
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        co2="data/eea/UNFCCC_v23.csv",
        swiss="data/switzerland-sfoe/switzerland-new_format.csv",
        idees="data/jrc-idees-2015",
        district_heat_share="data/district_heat_share.csv",
        eurostat=input_eurostat,
    output:
        energy_name=resources("energy_totals.csv"),
        co2_name=resources("co2_totals.csv"),
        transport_name=resources("transport_data.csv"),
    threads: 16
    resources:
        mem_mb=10000,
    log:
        logs("build_energy_totals.log"),
    benchmark:
        benchmarks("build_energy_totals")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_energy_totals.py"


rule build_biomass_potentials:
    params:
        biomass=config_provider("biomass"),
    input:
        enspreso_biomass=HTTP.remote(
            "https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/ENSPRESO/ENSPRESO_BIOMASS.xlsx",
            keep_local=True,
        ),
        nuts2="data/nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson",  # https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/#nuts21
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        nuts3_population=ancient("data/bundle/nama_10r_3popgdp.tsv.gz"),
        swiss_cantons=ancient("data/bundle/ch_cantons.csv"),
        swiss_population=ancient("data/bundle/je-e-21.03.02.xls"),
        country_shapes=resources("country_shapes.geojson"),
    output:
        biomass_potentials_all=resources(
            "biomass_potentials_all_s{simpl}_{clusters}.csv"
        ),
        biomass_potentials=resources("biomass_potentials_s{simpl}_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_biomass_potentials_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_biomass_potentials_s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_biomass_potentials.py"


if config["sector"]["biomass_transport"] or config["sector"]["biomass_spatial"]:

    rule build_biomass_transport_costs:
        input:
            transport_cost_data=HTTP.remote(
                "publications.jrc.ec.europa.eu/repository/bitstream/JRC98626/biomass potentials in europe_web rev.pdf",
                keep_local=True,
            ),
        output:
            biomass_transport_costs=resources("biomass_transport_costs.csv"),
        threads: 1
        resources:
            mem_mb=1000,
        log:
            logs("build_biomass_transport_costs.log"),
        benchmark:
            benchmarks("build_biomass_transport_costs")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_biomass_transport_costs.py"

    build_biomass_transport_costs_output = rules.build_biomass_transport_costs.output


if not (config["sector"]["biomass_transport"] or config["sector"]["biomass_spatial"]):
    # this is effecively an `else` statement which is however not liked by snakefmt
    build_biomass_transport_costs_output = {}


if config["sector"]["regional_co2_sequestration_potential"]["enable"]:

    rule build_sequestration_potentials:
        params:
            sequestration_potential=config_provider(
                "sector", "regional_co2_sequestration_potential"
            ),
        input:
            sequestration_potential=HTTP.remote(
                "https://raw.githubusercontent.com/ericzhou571/Co2Storage/main/resources/complete_map_2020_unit_Mt.geojson",
                keep_local=True,
            ),
            regions_onshore=resources(
                "regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
            regions_offshore=resources(
                "regions_offshore_elec_s{simpl}_{clusters}.geojson"
            ),
        output:
            sequestration_potential=resources(
                "co2_sequestration_potential_elec_s{simpl}_{clusters}.csv"
            ),
        threads: 1
        resources:
            mem_mb=4000,
        log:
            logs("build_sequestration_potentials_s{simpl}_{clusters}.log"),
        benchmark:
            benchmarks("build_sequestration_potentials_s{simpl}_{clusters}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_sequestration_potentials.py"

    build_sequestration_potentials_output = rules.build_sequestration_potentials.output


if not config["sector"]["regional_co2_sequestration_potential"]["enable"]:
    # this is effecively an `else` statement which is however not liked by snakefmt
    build_sequestration_potentials_output = {}


rule build_salt_cavern_potentials:
    input:
        salt_caverns="data/h2_salt_caverns_GWh_per_sqkm.geojson",
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
    output:
        h2_cavern_potential=resources("salt_cavern_potentials_s{simpl}_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_salt_cavern_potentials_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_salt_cavern_potentials_s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_salt_cavern_potentials.py"


rule build_ammonia_production:
    params:
        countries=config_provider("countries"),
    input:
        usgs="data/myb1-2017-nitro.xls",
    output:
        ammonia_production=resources("ammonia_production.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_ammonia_production.log"),
    benchmark:
        benchmarks("build_ammonia_production")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_ammonia_production.py"


rule build_industry_sector_ratios:
    params:
        industry=config_provider("industry"),
        ammonia=config_provider("sector", "ammonia", default=False),
    input:
        ammonia_production=resources("ammonia_production.csv"),
        idees="data/jrc-idees-2015",
    output:
        industry_sector_ratios=resources("industry_sector_ratios.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industry_sector_ratios.log"),
    benchmark:
        benchmarks("build_industry_sector_ratios")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industry_sector_ratios.py"


rule build_industrial_production_per_country:
    params:
        industry=config_provider("industry"),
        countries=config_provider("countries"),
    input:
        ammonia_production=resources("ammonia_production.csv"),
        jrc="data/jrc-idees-2015",
        eurostat="data/eurostat-energy_balances-may_2018_edition",
    output:
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    threads: 8
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_production_per_country.log"),
    benchmark:
        benchmarks("build_industrial_production_per_country")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_country.py"


rule build_industrial_production_per_country_tomorrow:
    params:
        industry=config_provider("industry"),
    input:
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    output:
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_production_per_country_tomorrow_{planning_horizons}.log"),
    benchmark:
        (
            benchmarks(
                "build_industrial_production_per_country_tomorrow_{planning_horizons}"
            )
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_country_tomorrow.py"


rule build_industrial_distribution_key:
    params:
        hotmaps_locate_missing=config_provider(
            "industry", "hotmaps_locate_missing", default=False
        ),
        countries=config_provider("countries"),
    input:
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
        hotmaps_industrial_database="data/Industrial_Database.csv",
    output:
        industrial_distribution_key=resources(
            "industrial_distribution_key_elec_s{simpl}_{clusters}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_distribution_key_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_industrial_distribution_key/s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_distribution_key.py"


rule build_industrial_production_per_node:
    input:
        industrial_distribution_key=resources(
            "industrial_distribution_key_elec_s{simpl}_{clusters}.csv"
        ),
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}.csv"
        ),
    output:
        industrial_production_per_node=resources(
            "industrial_production_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs(
            "build_industrial_production_per_node_s{simpl}_{clusters}_{planning_horizons}.log"
        ),
    benchmark:
        (
            benchmarks(
                "build_industrial_production_per_node/s{simpl}_{clusters}_{planning_horizons}"
            )
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_node.py"


rule build_industrial_energy_demand_per_node:
    input:
        industry_sector_ratios=resources("industry_sector_ratios.csv"),
        industrial_production_per_node=resources(
            "industrial_production_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
        industrial_energy_demand_per_node_today=resources(
            "industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv"
        ),
    output:
        industrial_energy_demand_per_node=resources(
            "industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs(
            "build_industrial_energy_demand_per_node_s{simpl}_{clusters}_{planning_horizons}.log"
        ),
    benchmark:
        (
            benchmarks(
                "build_industrial_energy_demand_per_node/s{simpl}_{clusters}_{planning_horizons}"
            )
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_node.py"


rule build_industrial_energy_demand_per_country_today:
    params:
        countries=config_provider("countries"),
        industry=config_provider("industry"),
    input:
        jrc="data/jrc-idees-2015",
        ammonia_production=resources("ammonia_production.csv"),
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    output:
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
    threads: 8
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_energy_demand_per_country_today.log"),
    benchmark:
        benchmarks("build_industrial_energy_demand_per_country_today")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_country_today.py"


rule build_industrial_energy_demand_per_node_today:
    input:
        industrial_distribution_key=resources(
            "industrial_distribution_key_elec_s{simpl}_{clusters}.csv"
        ),
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
    output:
        industrial_energy_demand_per_node_today=resources(
            "industrial_energy_demand_today_elec_s{simpl}_{clusters}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_energy_demand_per_node_today_s{simpl}_{clusters}.log"),
    benchmark:
        benchmarks("build_industrial_energy_demand_per_node_today/s{simpl}_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_node_today.py"


if config["sector"]["retrofitting"]["retro_endogen"]:

    rule build_retro_cost:
        params:
            retrofitting=config_provider("sector", "retrofitting"),
            countries=config_provider("countries"),
        input:
            building_stock="data/retro/data_building_stock.csv",
            data_tabula="data/retro/tabula-calculator-calcsetbuilding.csv",
            air_temperature=resources("temp_air_total_elec_s{simpl}_{clusters}.nc"),
            u_values_PL="data/retro/u_values_poland.csv",
            tax_w="data/retro/electricity_taxes_eu.csv",
            construction_index="data/retro/comparative_level_investment.csv",
            floor_area_missing="data/retro/floor_area_missing.csv",
            clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
            cost_germany="data/retro/retro_cost_germany.csv",
            window_assumptions="data/retro/window_assumptions.csv",
        output:
            retro_cost=resources("retro_cost_elec_s{simpl}_{clusters}.csv"),
            floor_area=resources("floor_area_elec_s{simpl}_{clusters}.csv"),
        resources:
            mem_mb=1000,
        log:
            logs("build_retro_cost_s{simpl}_{clusters}.log"),
        benchmark:
            benchmarks("build_retro_cost/s{simpl}_{clusters}")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_retro_cost.py"

    build_retro_cost_output = rules.build_retro_cost.output


if not config["sector"]["retrofitting"]["retro_endogen"]:
    # this is effecively an `else` statement which is however not liked by snakefmt
    build_retro_cost_output = {}


rule build_population_weighted_energy_totals:
    input:
        energy_totals=resources("energy_totals.csv"),
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
    output:
        resources("pop_weighted_energy_totals_s{simpl}_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_population_weighted_energy_totals_s{simpl}_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_population_weighted_energy_totals.py"


rule build_shipping_demand:
    input:
        ports="data/attributed_ports.json",
        scope=resources("europe_shape.geojson"),
        regions=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        demand=resources("energy_totals.csv"),
    output:
        resources("shipping_demand_s{simpl}_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_shipping_demand_s{simpl}_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_shipping_demand.py"


rule build_transport_demand:
    params:
        snapshots=config_provider("snapshots"),
        sector=config_provider("sector"),
    input:
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s{simpl}_{clusters}.csv"
        ),
        transport_data=resources("transport_data.csv"),
        traffic_data_KFZ="data/emobility/KFZ__count",
        traffic_data_Pkw="data/emobility/Pkw__count",
        temp_air_total=resources("temp_air_total_elec_s{simpl}_{clusters}.nc"),
    output:
        transport_demand=resources("transport_demand_s{simpl}_{clusters}.csv"),
        transport_data=resources("transport_data_s{simpl}_{clusters}.csv"),
        avail_profile=resources("avail_profile_s{simpl}_{clusters}.csv"),
        dsm_profile=resources("dsm_profile_s{simpl}_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_transport_demand_s{simpl}_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_transport_demand.py"


rule prepare_sector_network:
    params:
        co2_budget=config_provider("co2_budget"),
        conventional_carriers=config_provider(
            "existing_capacities", "conventional_carriers"
        ),
        foresight=config_provider("foresight"),
        costs=config_provider("costs"),
        sector=config_provider("sector"),
        industry=config_provider("industry"),
        pypsa_eur=config_provider("pypsa_eur"),
        length_factor=config_provider("lines", "length_factor"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        countries=config_provider("countries"),
        emissions_scope=config_provider("energy", "emissions"),
        eurostat_report_year=config_provider("energy", "eurostat_report_year"),
        RDIR=RDIR,
    input:
        **build_retro_cost_output,
        **build_biomass_transport_costs_output,
        **gas_infrastructure,
        **build_sequestration_potentials_output,
        network=resources("networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"),
        energy_totals_name=resources("energy_totals.csv"),
        eurostat=input_eurostat,
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s{simpl}_{clusters}.csv"
        ),
        shipping_demand=resources("shipping_demand_s{simpl}_{clusters}.csv"),
        transport_demand=resources("transport_demand_s{simpl}_{clusters}.csv"),
        transport_data=resources("transport_data_s{simpl}_{clusters}.csv"),
        avail_profile=resources("avail_profile_s{simpl}_{clusters}.csv"),
        dsm_profile=resources("dsm_profile_s{simpl}_{clusters}.csv"),
        co2_totals_name=resources("co2_totals.csv"),
        co2="data/eea/UNFCCC_v23.csv",
        biomass_potentials=resources("biomass_potentials_s{simpl}_{clusters}.csv"),
        heat_profile="data/heat_load_profile_BDEW.csv",
        costs="data/costs_{}.csv".format(config["costs"]["year"])
        if config["foresight"] == "overnight"
        else "data/costs_{planning_horizons}.csv",
        profile_offwind_ac=resources("profile_offwind-ac.nc"),
        profile_offwind_dc=resources("profile_offwind-dc.nc"),
        h2_cavern=resources("salt_cavern_potentials_s{simpl}_{clusters}.csv"),
        busmap_s=resources("busmap_elec_s{simpl}.csv"),
        busmap=resources("busmap_elec_s{simpl}_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
        simplified_pop_layout=resources("pop_layout_elec_s{simpl}.csv"),
        industrial_demand=resources(
            "industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
        heat_demand_urban=resources("heat_demand_urban_elec_s{simpl}_{clusters}.nc"),
        heat_demand_rural=resources("heat_demand_rural_elec_s{simpl}_{clusters}.nc"),
        heat_demand_total=resources("heat_demand_total_elec_s{simpl}_{clusters}.nc"),
        temp_soil_total=resources("temp_soil_total_elec_s{simpl}_{clusters}.nc"),
        temp_soil_rural=resources("temp_soil_rural_elec_s{simpl}_{clusters}.nc"),
        temp_soil_urban=resources("temp_soil_urban_elec_s{simpl}_{clusters}.nc"),
        temp_air_total=resources("temp_air_total_elec_s{simpl}_{clusters}.nc"),
        temp_air_rural=resources("temp_air_rural_elec_s{simpl}_{clusters}.nc"),
        temp_air_urban=resources("temp_air_urban_elec_s{simpl}_{clusters}.nc"),
        cop_soil_total=resources("cop_soil_total_elec_s{simpl}_{clusters}.nc"),
        cop_soil_rural=resources("cop_soil_rural_elec_s{simpl}_{clusters}.nc"),
        cop_soil_urban=resources("cop_soil_urban_elec_s{simpl}_{clusters}.nc"),
        cop_air_total=resources("cop_air_total_elec_s{simpl}_{clusters}.nc"),
        cop_air_rural=resources("cop_air_rural_elec_s{simpl}_{clusters}.nc"),
        cop_air_urban=resources("cop_air_urban_elec_s{simpl}_{clusters}.nc"),
        solar_thermal_total=resources(
            "solar_thermal_total_elec_s{simpl}_{clusters}.nc"
        )
        if config["sector"]["solar_thermal"]
        else [],
        solar_thermal_urban=resources(
            "solar_thermal_urban_elec_s{simpl}_{clusters}.nc"
        )
        if config["sector"]["solar_thermal"]
        else [],
        solar_thermal_rural=resources(
            "solar_thermal_rural_elec_s{simpl}_{clusters}.nc"
        )
        if config["sector"]["solar_thermal"]
        else [],
    output:
        RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS
        + "prepare_sector_network_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "prepare_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/prepare_sector_network.py"
