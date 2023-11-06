# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule build_population_layouts:
    input:
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
        urban_percent="data/urban_percent.csv",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        pop_layout_total=RESOURCES + "pop_layout_total.nc",
        pop_layout_urban=RESOURCES + "pop_layout_urban.nc",
        pop_layout_rural=RESOURCES + "pop_layout_rural.nc",
    log:
        LOGS + "build_population_layouts.log",
    resources:
        mem_mb=20000,
    benchmark:
        BENCHMARKS + "build_population_layouts"
    threads: 8
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total=RESOURCES + "pop_layout_total.nc",
        pop_layout_urban=RESOURCES + "pop_layout_urban.nc",
        pop_layout_rural=RESOURCES + "pop_layout_rural.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
    log:
        LOGS + "build_clustered_population_layouts_{simpl}_{clusters}.log",
    resources:
        mem_mb=10000,
    benchmark:
        BENCHMARKS + "build_clustered_population_layouts/s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_clustered_population_layouts.py"


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
    log:
        LOGS + "build_simplified_population_layouts_{simpl}",
    benchmark:
        BENCHMARKS + "build_simplified_population_layouts/s{simpl}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_clustered_population_layouts.py"


if config["sector"]["gas_network"] or config["sector"]["H2_retrofit"]:

    rule build_gas_network:
        input:
            gas_network="data/gas_network/scigrid-gas/data/IGGIELGN_PipeSegments.geojson",
        output:
            cleaned_gas_network=RESOURCES + "gas_network.csv",
        resources:
            mem_mb=4000,
        log:
            LOGS + "build_gas_network.log",
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
        log:
            LOGS + "build_gas_input_locations_s{simpl}_{clusters}.log",
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_gas_input_locations.py"

    rule cluster_gas_network:
        input:
            cleaned_gas_network=RESOURCES + "gas_network.csv",
            regions_onshore=RESOURCES
            + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
            regions_offshore=RESOURCES
            + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
        output:
            clustered_gas_network=RESOURCES + "gas_network_elec_s{simpl}_{clusters}.csv",
        resources:
            mem_mb=4000,
        log:
            LOGS + "cluster_gas_network_s{simpl}_{clusters}.log",
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
        snapshots=config["snapshots"],
    input:
        pop_layout=RESOURCES + "pop_layout_{scope}.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        heat_demand=RESOURCES + "heat_demand_{scope}_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    threads: 8
    log:
        LOGS + "build_heat_demands_{scope}_{simpl}_{clusters}.loc",
    benchmark:
        BENCHMARKS + "build_heat_demands/{scope}_s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_heat_demand.py"


rule build_temperature_profiles:
    params:
        snapshots=config["snapshots"],
    input:
        pop_layout=RESOURCES + "pop_layout_{scope}.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        temp_soil=RESOURCES + "temp_soil_{scope}_elec_s{simpl}_{clusters}.nc",
        temp_air=RESOURCES + "temp_air_{scope}_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    threads: 8
    log:
        LOGS + "build_temperature_profiles_{scope}_{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "build_temperature_profiles/{scope}_s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_temperature_profiles.py"


rule build_cop_profiles:
    params:
        heat_pump_sink_T=config["sector"]["heat_pump_sink_T"],
    input:
        temp_soil_total=RESOURCES + "temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural=RESOURCES + "temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban=RESOURCES + "temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total=RESOURCES + "temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural=RESOURCES + "temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban=RESOURCES + "temp_air_urban_elec_s{simpl}_{clusters}.nc",
    output:
        cop_soil_total=RESOURCES + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural=RESOURCES + "cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban=RESOURCES + "cop_soil_urban_elec_s{simpl}_{clusters}.nc",
        cop_air_total=RESOURCES + "cop_air_total_elec_s{simpl}_{clusters}.nc",
        cop_air_rural=RESOURCES + "cop_air_rural_elec_s{simpl}_{clusters}.nc",
        cop_air_urban=RESOURCES + "cop_air_urban_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    log:
        LOGS + "build_cop_profiles_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "build_cop_profiles/s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_cop_profiles.py"


rule build_solar_thermal_profiles:
    params:
        snapshots=config["snapshots"],
        solar_thermal=config["solar_thermal"],
    input:
        pop_layout=RESOURCES + "pop_layout_{scope}.nc",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        cutout="cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
    output:
        solar_thermal=RESOURCES + "solar_thermal_{scope}_elec_s{simpl}_{clusters}.nc",
    resources:
        mem_mb=20000,
    threads: 16
    log:
        LOGS + "build_solar_thermal_profiles_{scope}_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "build_solar_thermal_profiles/{scope}_s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_solar_thermal_profiles.py"


rule build_energy_totals:
    params:
        countries=config["countries"],
        energy=config["energy"],
    input:
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
        co2="data/bundle-sector/eea/UNFCCC_v23.csv",
        swiss="data/bundle-sector/switzerland-sfoe/switzerland-new_format.csv",
        idees="data/bundle-sector/jrc-idees-2015",
        district_heat_share="data/district_heat_share.csv",
        eurostat=input_eurostat,
    output:
        energy_name=RESOURCES + "energy_totals.csv",
        co2_name=RESOURCES + "co2_totals.csv",
        transport_name=RESOURCES + "transport_data.csv",
    threads: 16
    resources:
        mem_mb=10000,
    log:
        LOGS + "build_energy_totals.log",
    benchmark:
        BENCHMARKS + "build_energy_totals"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_energy_totals.py"


rule build_biomass_potentials:
    params:
        biomass=config["biomass"],
    input:
        enspreso_biomass=HTTP.remote(
            "https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/ENSPRESO/ENSPRESO_BIOMASS.xlsx",
            keep_local=True,
        ),
        nuts2="data/bundle-sector/nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson",  # https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/#nuts21
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        nuts3_population=ancient("data/bundle/nama_10r_3popgdp.tsv.gz"),
        swiss_cantons=ancient("data/bundle/ch_cantons.csv"),
        swiss_population=ancient("data/bundle/je-e-21.03.02.xls"),
        country_shapes=RESOURCES + "country_shapes.geojson",
    output:
        biomass_potentials_all=RESOURCES
        + "biomass_potentials_all_s{simpl}_{clusters}_{planning_horizons}.csv",
        biomass_potentials=RESOURCES
        + "biomass_potentials_s{simpl}_{clusters}_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        LOGS + "build_biomass_potentials_s{simpl}_{clusters}_{planning_horizons}.log",
    benchmark:
        BENCHMARKS + "build_biomass_potentials_s{simpl}_{clusters}_{planning_horizons}"
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
            biomass_transport_costs=RESOURCES + "biomass_transport_costs.csv",
        threads: 1
        resources:
            mem_mb=1000,
        log:
            LOGS + "build_biomass_transport_costs.log",
        benchmark:
            BENCHMARKS + "build_biomass_transport_costs"
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
            sequestration_potential=config["sector"][
                "regional_co2_sequestration_potential"
            ],
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
        log:
            LOGS + "build_sequestration_potentials_s{simpl}_{clusters}.log",
        benchmark:
            BENCHMARKS + "build_sequestration_potentials_s{simpl}_{clusters}"
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
        salt_caverns="data/bundle-sector/h2_salt_caverns_GWh_per_sqkm.geojson",
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        regions_offshore=RESOURCES + "regions_offshore_elec_s{simpl}_{clusters}.geojson",
    output:
        h2_cavern_potential=RESOURCES + "salt_cavern_potentials_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS + "build_salt_cavern_potentials_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "build_salt_cavern_potentials_s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_salt_cavern_potentials.py"


rule build_ammonia_production:
    params:
        countries=config["countries"],
    input:
        usgs="data/bundle-sector/myb1-2017-nitro.xls",
    output:
        ammonia_production=RESOURCES + "ammonia_production.csv",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        LOGS + "build_ammonia_production.log",
    benchmark:
        BENCHMARKS + "build_ammonia_production"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_ammonia_production.py"


rule build_industry_sector_ratios:
    params:
        industry=config["industry"],
        ammonia=config["sector"].get("ammonia", False),
    input:
        ammonia_production=RESOURCES + "ammonia_production.csv",
        idees="data/bundle-sector/jrc-idees-2015",
    output:
        industry_sector_ratios=RESOURCES + "industry_sector_ratios.csv",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        LOGS + "build_industry_sector_ratios.log",
    benchmark:
        BENCHMARKS + "build_industry_sector_ratios"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industry_sector_ratios.py"


rule build_industrial_production_per_country:
    params:
        industry=config["industry"],
        countries=config["countries"],
    input:
        ammonia_production=RESOURCES + "ammonia_production.csv",
        jrc="data/bundle-sector/jrc-idees-2015",
        eurostat="data/bundle-sector/eurostat-energy_balances-may_2018_edition",
    output:
        industrial_production_per_country=RESOURCES
        + "industrial_production_per_country.csv",
    threads: 8
    resources:
        mem_mb=1000,
    log:
        LOGS + "build_industrial_production_per_country.log",
    benchmark:
        BENCHMARKS + "build_industrial_production_per_country"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_country.py"


rule build_industrial_production_per_country_tomorrow:
    params:
        industry=config["industry"],
    input:
        industrial_production_per_country=RESOURCES
        + "industrial_production_per_country.csv",
    output:
        industrial_production_per_country_tomorrow=RESOURCES
        + "industrial_production_per_country_tomorrow_{planning_horizons}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        LOGS
        + "build_industrial_production_per_country_tomorrow_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "build_industrial_production_per_country_tomorrow_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_country_tomorrow.py"


rule build_industrial_distribution_key:
    params:
        hotmaps_locate_missing=config["industry"].get("hotmaps_locate_missing", False),
        countries=config["countries"],
    input:
        regions_onshore=RESOURCES + "regions_onshore_elec_s{simpl}_{clusters}.geojson",
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
        hotmaps_industrial_database="data/bundle-sector/Industrial_Database.csv",
    output:
        industrial_distribution_key=RESOURCES
        + "industrial_distribution_key_elec_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        LOGS + "build_industrial_distribution_key_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "build_industrial_distribution_key/s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_distribution_key.py"


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
    log:
        LOGS
        + "build_industrial_production_per_node_s{simpl}_{clusters}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "build_industrial_production_per_node/s{simpl}_{clusters}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_node.py"


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
    log:
        LOGS
        + "build_industrial_energy_demand_per_node_s{simpl}_{clusters}_{planning_horizons}.log",
    benchmark:
        (
            BENCHMARKS
            + "build_industrial_energy_demand_per_node/s{simpl}_{clusters}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_node.py"


rule build_industrial_energy_demand_per_country_today:
    params:
        countries=config["countries"],
        industry=config["industry"],
    input:
        jrc="data/bundle-sector/jrc-idees-2015",
        ammonia_production=RESOURCES + "ammonia_production.csv",
        industrial_production_per_country=RESOURCES
        + "industrial_production_per_country.csv",
    output:
        industrial_energy_demand_per_country_today=RESOURCES
        + "industrial_energy_demand_per_country_today.csv",
    threads: 8
    resources:
        mem_mb=1000,
    log:
        LOGS + "build_industrial_energy_demand_per_country_today.log",
    benchmark:
        BENCHMARKS + "build_industrial_energy_demand_per_country_today"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_country_today.py"


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
    log:
        LOGS + "build_industrial_energy_demand_per_node_today_s{simpl}_{clusters}.log",
    benchmark:
        BENCHMARKS + "build_industrial_energy_demand_per_node_today/s{simpl}_{clusters}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_node_today.py"


if config["sector"]["retrofitting"]["retro_endogen"]:

    rule build_retro_cost:
        params:
            retrofitting=config["sector"]["retrofitting"],
            countries=config["countries"],
        input:
            building_stock="data/retro/data_building_stock.csv",
            data_tabula="data/bundle-sector/retro/tabula-calculator-calcsetbuilding.csv",
            air_temperature=RESOURCES + "temp_air_total_elec_s{simpl}_{clusters}.nc",
            u_values_PL="data/retro/u_values_poland.csv",
            tax_w="data/retro/electricity_taxes_eu.csv",
            construction_index="data/retro/comparative_level_investment.csv",
            floor_area_missing="data/retro/floor_area_missing.csv",
            clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
            cost_germany="data/retro/retro_cost_germany.csv",
            window_assumptions="data/retro/window_assumptions.csv",
        output:
            retro_cost=RESOURCES + "retro_cost_elec_s{simpl}_{clusters}.csv",
            floor_area=RESOURCES + "floor_area_elec_s{simpl}_{clusters}.csv",
        resources:
            mem_mb=1000,
        log:
            LOGS + "build_retro_cost_s{simpl}_{clusters}.log",
        benchmark:
            BENCHMARKS + "build_retro_cost/s{simpl}_{clusters}"
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
        energy_totals=RESOURCES + "energy_totals.csv",
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
    output:
        RESOURCES + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS + "build_population_weighted_energy_totals_s{simpl}_{clusters}.log",
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_population_weighted_energy_totals.py"


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
    log:
        LOGS + "build_shipping_demand_s{simpl}_{clusters}.log",
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_shipping_demand.py"


rule build_transport_demand:
    params:
        snapshots=config["snapshots"],
        sector=config["sector"],
    input:
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
        pop_weighted_energy_totals=RESOURCES
        + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        transport_data=RESOURCES + "transport_data.csv",
        traffic_data_KFZ="data/bundle-sector/emobility/KFZ__count",
        traffic_data_Pkw="data/bundle-sector/emobility/Pkw__count",
        temp_air_total=RESOURCES + "temp_air_total_elec_s{simpl}_{clusters}.nc",
    output:
        transport_demand=RESOURCES + "transport_demand_s{simpl}_{clusters}.csv",
        transport_data=RESOURCES + "transport_data_s{simpl}_{clusters}.csv",
        avail_profile=RESOURCES + "avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile=RESOURCES + "dsm_profile_s{simpl}_{clusters}.csv",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        LOGS + "build_transport_demand_s{simpl}_{clusters}.log",
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_transport_demand.py"


rule prepare_sector_network:
    params:
        co2_budget=config["co2_budget"],
        conventional_carriers=config["existing_capacities"]["conventional_carriers"],
        foresight=config["foresight"],
        costs=config["costs"],
        sector=config["sector"],
        industry=config["industry"],
        pypsa_eur=config["pypsa_eur"],
        length_factor=config["lines"]["length_factor"],
        planning_horizons=config["scenario"]["planning_horizons"],
        countries=config["countries"],
        emissions_scope=config["energy"]["emissions"],
        eurostat_report_year=config["energy"]["eurostat_report_year"],
        RDIR=RDIR,
    input:
        **build_retro_cost_output,
        **build_biomass_transport_costs_output,
        **gas_infrastructure,
        **build_sequestration_potentials_output,
        network=RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        energy_totals_name=RESOURCES + "energy_totals.csv",
        eurostat=input_eurostat,
        pop_weighted_energy_totals=RESOURCES
        + "pop_weighted_energy_totals_s{simpl}_{clusters}.csv",
        shipping_demand=RESOURCES + "shipping_demand_s{simpl}_{clusters}.csv",
        transport_demand=RESOURCES + "transport_demand_s{simpl}_{clusters}.csv",
        transport_data=RESOURCES + "transport_data_s{simpl}_{clusters}.csv",
        avail_profile=RESOURCES + "avail_profile_s{simpl}_{clusters}.csv",
        dsm_profile=RESOURCES + "dsm_profile_s{simpl}_{clusters}.csv",
        co2_totals_name=RESOURCES + "co2_totals.csv",
        co2="data/bundle-sector/eea/UNFCCC_v23.csv",
        biomass_potentials=RESOURCES
        + "biomass_potentials_s{simpl}_{clusters}_"
        + "{}.csv".format(config["biomass"]["year"])
        if config["foresight"] == "overnight"
        else RESOURCES
        + "biomass_potentials_s{simpl}_{clusters}_{planning_horizons}.csv",
        heat_profile="data/heat_load_profile_BDEW.csv",
        costs="data/costs_{}.csv".format(config["costs"]["year"])
        if config["foresight"] == "overnight"
        else "data/costs_{planning_horizons}.csv",
        profile_offwind_ac=RESOURCES + "profile_offwind-ac.nc",
        profile_offwind_dc=RESOURCES + "profile_offwind-dc.nc",
        profile_offwind_float=RESOURCES + "profile_offwind-float.nc",
        h2_cavern=RESOURCES + "salt_cavern_potentials_s{simpl}_{clusters}.csv",
        busmap_s=RESOURCES + "busmap_elec_s{simpl}.csv",
        busmap=RESOURCES + "busmap_elec_s{simpl}_{clusters}.csv",
        clustered_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}_{clusters}.csv",
        simplified_pop_layout=RESOURCES + "pop_layout_elec_s{simpl}.csv",
        industrial_demand=RESOURCES
        + "industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv",
        heat_demand_urban=RESOURCES + "heat_demand_urban_elec_s{simpl}_{clusters}.nc",
        heat_demand_rural=RESOURCES + "heat_demand_rural_elec_s{simpl}_{clusters}.nc",
        heat_demand_total=RESOURCES + "heat_demand_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_total=RESOURCES + "temp_soil_total_elec_s{simpl}_{clusters}.nc",
        temp_soil_rural=RESOURCES + "temp_soil_rural_elec_s{simpl}_{clusters}.nc",
        temp_soil_urban=RESOURCES + "temp_soil_urban_elec_s{simpl}_{clusters}.nc",
        temp_air_total=RESOURCES + "temp_air_total_elec_s{simpl}_{clusters}.nc",
        temp_air_rural=RESOURCES + "temp_air_rural_elec_s{simpl}_{clusters}.nc",
        temp_air_urban=RESOURCES + "temp_air_urban_elec_s{simpl}_{clusters}.nc",
        cop_soil_total=RESOURCES + "cop_soil_total_elec_s{simpl}_{clusters}.nc",
        cop_soil_rural=RESOURCES + "cop_soil_rural_elec_s{simpl}_{clusters}.nc",
        cop_soil_urban=RESOURCES + "cop_soil_urban_elec_s{simpl}_{clusters}.nc",
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
