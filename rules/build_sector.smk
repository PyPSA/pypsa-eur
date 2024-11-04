# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule build_population_layouts:
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        urban_percent="data/worldbank/API_SP.URB.TOTL.IN.ZS_DS2_en_csv_v2.csv",
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("atlite", "default_cutout")(w)
        + ".nc",
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
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("atlite", "default_cutout")(w)
        + ".nc",
    output:
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
    log:
        logs("build_clustered_population_layouts_s_{clusters}.log"),
    resources:
        mem_mb=10000,
    benchmark:
        benchmarks("build_clustered_population_layouts/s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_clustered_population_layouts.py"


rule build_simplified_population_layouts:
    input:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
        regions_onshore=resources("regions_onshore_base_s.geojson"),
        cutout=lambda w: "cutouts/"
        + CDIR
        + config_provider("atlite", "default_cutout")(w)
        + ".nc",
    output:
        clustered_pop_layout=resources("pop_layout_base_s.csv"),
    resources:
        mem_mb=10000,
    log:
        logs("build_simplified_population_layouts_s"),
    benchmark:
        benchmarks("build_simplified_population_layouts/s")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_clustered_population_layouts.py"


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
        gem="data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        entry="data/gas_network/scigrid-gas/data/IGGIELGN_BorderPoints.geojson",
        storage="data/gas_network/scigrid-gas/data/IGGIELGN_Storages.geojson",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
    output:
        gas_input_nodes=resources("gas_input_locations_s_{clusters}.geojson"),
        gas_input_nodes_simplified=resources(
            "gas_input_locations_s_{clusters}_simplified.csv"
        ),
    resources:
        mem_mb=2000,
    log:
        logs("build_gas_input_locations_s_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_gas_input_locations.py"


rule cluster_gas_network:
    input:
        cleaned_gas_network=resources("gas_network.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
    output:
        clustered_gas_network=resources("gas_network_base_s_{clusters}.csv"),
    resources:
        mem_mb=4000,
    log:
        logs("cluster_gas_network_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/cluster_gas_network.py"


def heat_demand_cutout(wildcards):
    c = config_provider("sector", "heat_demand_cutout")(wildcards)
    if c == "default":
        return (
            "cutouts/"
            + CDIR
            + config_provider("atlite", "default_cutout")(wildcards)
            + ".nc"
        )
    else:
        return "cutouts/" + CDIR + c + ".nc"


rule build_daily_heat_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        pop_layout=resources("pop_layout_total.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=heat_demand_cutout,
    output:
        heat_demand=resources("daily_heat_demand_total_base_s_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 8
    log:
        logs("build_daily_heat_demand_total_s_{clusters}.loc"),
    benchmark:
        benchmarks("build_daily_heat_demand/total_s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_daily_heat_demand.py"


rule build_hourly_heat_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        heat_profile="data/heat_load_profile_BDEW.csv",
        heat_demand=resources("daily_heat_demand_total_base_s_{clusters}.nc"),
    output:
        heat_demand=resources("hourly_heat_demand_total_base_s_{clusters}.nc"),
    resources:
        mem_mb=2000,
    threads: 8
    log:
        logs("build_hourly_heat_demand_total_s_{clusters}.loc"),
    benchmark:
        benchmarks("build_hourly_heat_demand/total_s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_hourly_heat_demand.py"


rule build_temperature_profiles:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        pop_layout=resources("pop_layout_total.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=heat_demand_cutout,
    output:
        temp_soil=resources("temp_soil_total_base_s_{clusters}.nc"),
        temp_air=resources("temp_air_total_base_s_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 8
    log:
        logs("build_temperature_profiles_total_s_{clusters}.log"),
    benchmark:
        benchmarks("build_temperature_profiles/total_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_temperature_profiles.py"


rule build_central_heating_temperature_profiles:
    params:
        max_forward_temperature_central_heating_baseyear=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "max_forward_temperature_baseyear",
        ),
        min_forward_temperature_central_heating_baseyear=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "min_forward_temperature_baseyear",
        ),
        return_temperature_central_heating_baseyear=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "return_temperature_baseyear",
        ),
        snapshots=config_provider("snapshots"),
        lower_threshold_ambient_temperature=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "lower_threshold_ambient_temperature",
        ),
        upper_threshold_ambient_temperature=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "upper_threshold_ambient_temperature",
        ),
        rolling_window_ambient_temperature=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "rolling_window_ambient_temperature",
        ),
        relative_annual_temperature_reduction=config_provider(
            "sector",
            "district_heating",
            "supply_temperature_approximation",
            "relative_annual_temperature_reduction",
        ),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        temp_air_total=resources("temp_air_total_base_s_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
    resources:
        mem_mb=20000,
    log:
        logs(
            "build_central_heating_temperature_profiles_s_{clusters}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "build_central_heating_temperature_profiles/s_{clusters}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_central_heating_temperature_profiles/run.py"


rule build_cop_profiles:
    params:
        heat_pump_sink_T_decentral_heating=config_provider(
            "sector", "heat_pump_sink_T_individual_heating"
        ),
        heat_source_cooling_central_heating=config_provider(
            "sector", "district_heating", "heat_source_cooling"
        ),
        heat_pump_cop_approximation_central_heating=config_provider(
            "sector", "district_heating", "heat_pump_cop_approximation"
        ),
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        snapshots=config_provider("snapshots"),
    input:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        temp_soil_total=resources("temp_soil_total_base_s_{clusters}.nc"),
        temp_air_total=resources("temp_air_total_base_s_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
    resources:
        mem_mb=20000,
    log:
        logs("build_cop_profiles_s_{clusters}_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_cop_profiles/s_{clusters}_{planning_horizons}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_cop_profiles/run.py"


def solar_thermal_cutout(wildcards):
    c = config_provider("solar_thermal", "cutout")(wildcards)
    if c == "default":
        return (
            "cutouts/"
            + CDIR
            + config_provider("atlite", "default_cutout")(wildcards)
            + ".nc"
        )
    else:
        return "cutouts/" + CDIR + c + ".nc"


rule build_solar_thermal_profiles:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        solar_thermal=config_provider("solar_thermal"),
    input:
        pop_layout=resources("pop_layout_total.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=solar_thermal_cutout,
    output:
        solar_thermal=resources("solar_thermal_total_base_s_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 16
    log:
        logs("build_solar_thermal_profiles_total_s_{clusters}.log"),
    benchmark:
        benchmarks("build_solar_thermal_profiles/total_{clusters}")
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
        co2="data/bundle/eea/UNFCCC_v23.csv",
        swiss="data/switzerland-new_format-all_years.csv",
        swiss_transport="data/gr-e-11.03.02.01.01-cc.csv",
        idees="data/jrc-idees-2021",
        district_heat_share="data/district_heat_share.csv",
        eurostat="data/eurostat/Balances-April2023",
        eurostat_households="data/eurostat/eurostat-household_energy_balances-february_2024.csv",
    output:
        transformation_output_coke=resources("transformation_output_coke.csv"),
        energy_name=resources("energy_totals.csv"),
        co2_name=resources("co2_totals.csv"),
        transport_name=resources("transport_data.csv"),
        district_heat_share=resources("district_heat_share.csv"),
        heating_efficiencies=resources("heating_efficiencies.csv"),
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


rule build_heat_totals:
    input:
        hdd="data/era5-annual-HDD-per-country.csv",
        energy_totals=resources("energy_totals.csv"),
    output:
        heat_totals=resources("heat_totals.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_heat_totals.log"),
    benchmark:
        benchmarks("build_heat_totals")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_heat_totals.py"


rule build_biomass_potentials:
    params:
        biomass=config_provider("biomass"),
    input:
        enspreso_biomass="data/ENSPRESO_BIOMASS.xlsx",
        eurostat="data/eurostat/Balances-April2023",
        nuts2="data/nuts/NUTS_RG_03M_2013_4326_LEVL_2.geojson",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        nuts3_population=ancient("data/bundle/nama_10r_3popgdp.tsv.gz"),
        swiss_cantons=ancient("data/ch_cantons.csv"),
        swiss_population=ancient("data/bundle/je-e-21.03.02.xls"),
        country_shapes=resources("country_shapes.geojson"),
    output:
        biomass_potentials_all=resources(
            "biomass_potentials_all_{clusters}_{planning_horizons}.csv"
        ),
        biomass_potentials=resources(
            "biomass_potentials_s_{clusters}_{planning_horizons}.csv"
        ),
    threads: 8
    resources:
        mem_mb=1000,
    log:
        logs("build_biomass_potentials_s_{clusters}_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_biomass_potentials_s_{clusters}_{planning_horizons}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_biomass_potentials.py"


rule build_biomass_transport_costs:
    input:
        sc1="data/biomass_transport_costs_supplychain1.csv",
        sc2="data/biomass_transport_costs_supplychain2.csv",
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


rule build_sequestration_potentials:
    params:
        sequestration_potential=config_provider(
            "sector", "regional_co2_sequestration_potential"
        ),
    input:
        sequestration_potential="data/complete_map_2020_unit_Mt.geojson",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
    output:
        sequestration_potential=resources(
            "co2_sequestration_potential_base_s_{clusters}.csv"
        ),
    threads: 1
    resources:
        mem_mb=4000,
    log:
        logs("build_sequestration_potentials_{clusters}.log"),
    benchmark:
        benchmarks("build_sequestration_potentials_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_sequestration_potentials.py"


rule build_salt_cavern_potentials:
    input:
        salt_caverns="data/bundle/h2_salt_caverns_GWh_per_sqkm.geojson",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
    output:
        h2_cavern_potential=resources("salt_cavern_potentials_s_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_salt_cavern_potentials_s_{clusters}.log"),
    benchmark:
        benchmarks("build_salt_cavern_potentials_s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_salt_cavern_potentials.py"


rule build_ammonia_production:
    input:
        usgs="data/myb1-2022-nitro-ert.xlsx",
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
        idees="data/jrc-idees-2021",
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


rule build_industry_sector_ratios_intermediate:
    params:
        industry=config_provider("industry"),
    input:
        industry_sector_ratios=resources("industry_sector_ratios.csv"),
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    output:
        industry_sector_ratios=resources(
            "industry_sector_ratios_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industry_sector_ratios_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_industry_sector_ratios_{planning_horizons}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industry_sector_ratios_intermediate.py"


rule build_industrial_production_per_country:
    params:
        industry=config_provider("industry"),
        countries=config_provider("countries"),
    input:
        ch_industrial_production="data/ch_industrial_production_per_subsector.csv",
        ammonia_production=resources("ammonia_production.csv"),
        jrc="data/jrc-idees-2021",
        eurostat="data/eurostat/Balances-April2023",
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
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        hotmaps="data/Industrial_Database.csv",
        gem_gspt="data/gem/Global-Steel-Plant-Tracker-April-2024-Standard-Copy-V1.xlsx",
        ammonia="data/ammonia_plants.csv",
        cement_supplement="data/cement-plants-noneu.csv",
        refineries_supplement="data/refineries-noneu.csv",
    output:
        industrial_distribution_key=resources(
            "industrial_distribution_key_base_s_{clusters}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_distribution_key_{clusters}.log"),
    benchmark:
        benchmarks("build_industrial_distribution_key/s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_distribution_key.py"


rule build_industrial_production_per_node:
    input:
        industrial_distribution_key=resources(
            "industrial_distribution_key_base_s_{clusters}.csv"
        ),
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{planning_horizons}.csv"
        ),
    output:
        industrial_production_per_node=resources(
            "industrial_production_base_s_{clusters}_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_production_per_node_{clusters}_{planning_horizons}.log"),
    benchmark:
        (
            benchmarks(
                "build_industrial_production_per_node/s_{clusters}_{planning_horizons}"
            )
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_production_per_node.py"


rule build_industrial_energy_demand_per_node:
    input:
        industry_sector_ratios=resources(
            "industry_sector_ratios_{planning_horizons}.csv"
        ),
        industrial_production_per_node=resources(
            "industrial_production_base_s_{clusters}_{planning_horizons}.csv"
        ),
        industrial_energy_demand_per_node_today=resources(
            "industrial_energy_demand_today_base_s_{clusters}.csv"
        ),
    output:
        industrial_energy_demand_per_node=resources(
            "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs(
            "build_industrial_energy_demand_per_node_{clusters}_{planning_horizons}.log"
        ),
    benchmark:
        (
            benchmarks(
                "build_industrial_energy_demand_per_node/s_{clusters}_{planning_horizons}"
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
        ammonia=config_provider("sector", "ammonia", default=False),
    input:
        transformation_output_coke=resources("transformation_output_coke.csv"),
        jrc="data/jrc-idees-2021",
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
            "industrial_distribution_key_base_s_{clusters}.csv"
        ),
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
    output:
        industrial_energy_demand_per_node_today=resources(
            "industrial_energy_demand_today_base_s_{clusters}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industrial_energy_demand_per_node_today_{clusters}.log"),
    benchmark:
        benchmarks("build_industrial_energy_demand_per_node_today/s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_industrial_energy_demand_per_node_today.py"


rule build_retro_cost:
    params:
        retrofitting=config_provider("sector", "retrofitting"),
        countries=config_provider("countries"),
    input:
        building_stock="data/retro/data_building_stock.csv",
        data_tabula="data/bundle/retro/tabula-calculator-calcsetbuilding.csv",
        air_temperature=resources("temp_air_total_base_s_{clusters}.nc"),
        u_values_PL="data/retro/u_values_poland.csv",
        tax_w="data/retro/electricity_taxes_eu.csv",
        construction_index="data/retro/comparative_level_investment.csv",
        floor_area_missing="data/retro/floor_area_missing.csv",
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        cost_germany="data/retro/retro_cost_germany.csv",
        window_assumptions="data/retro/window_assumptions.csv",
    output:
        retro_cost=resources("retro_cost_base_s_{clusters}.csv"),
        floor_area=resources("floor_area_base_s_{clusters}.csv"),
    resources:
        mem_mb=1000,
    log:
        logs("build_retro_cost_{clusters}.log"),
    benchmark:
        benchmarks("build_retro_cost/s_{clusters}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_retro_cost.py"


rule build_population_weighted_energy_totals:
    params:
        snapshots=config_provider("snapshots"),
    input:
        energy_totals=resources("{kind}_totals.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
    output:
        resources("pop_weighted_{kind}_totals_s_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_population_weighted_{kind}_totals_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_population_weighted_energy_totals.py"


rule build_shipping_demand:
    input:
        ports="data/attributed_ports.json",
        scope=resources("europe_shape.geojson"),
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        demand=resources("energy_totals.csv"),
    params:
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    output:
        resources("shipping_demand_s_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_shipping_demand_s_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_shipping_demand.py"


rule build_transport_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        sector=config_provider("sector"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        transport_data=resources("transport_data.csv"),
        traffic_data_KFZ="data/bundle/emobility/KFZ__count",
        traffic_data_Pkw="data/bundle/emobility/Pkw__count",
        temp_air_total=resources("temp_air_total_base_s_{clusters}.nc"),
    output:
        transport_demand=resources("transport_demand_s_{clusters}.csv"),
        transport_data=resources("transport_data_s_{clusters}.csv"),
        avail_profile=resources("avail_profile_s_{clusters}.csv"),
        dsm_profile=resources("dsm_profile_s_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_transport_demand_s_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_transport_demand.py"


rule build_district_heat_share:
    params:
        sector=config_provider("sector"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        district_heat_share=resources("district_heat_share.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
    output:
        district_heat_share=resources(
            "district_heat_share_base_s_{clusters}_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_district_heat_share_{clusters}_{planning_horizons}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_district_heat_share.py"


rule build_existing_heating_distribution:
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
    input:
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        clustered_pop_energy_layout=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        district_heat_share=resources(
            "district_heat_share_base_s_{clusters}_{planning_horizons}.csv"
        ),
    output:
        existing_heating_distribution=resources(
            "existing_heating_distribution_base_s_{clusters}_{planning_horizons}.csv"
        ),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs(
            "build_existing_heating_distribution_base_s_{clusters}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "build_existing_heating_distribution/base_s_{clusters}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_existing_heating_distribution.py"


rule time_aggregation:
    params:
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        solver_name=config_provider("solving", "solver", "name"),
    input:
        network=resources("networks/base_s_{clusters}_elec_l{ll}_{opts}.nc"),
        hourly_heat_demand_total=lambda w: (
            resources("hourly_heat_demand_total_base_s_{clusters}.nc")
            if config_provider("sector", "heating")(w)
            else []
        ),
        solar_thermal_total=lambda w: (
            resources("solar_thermal_total_base_s_{clusters}.nc")
            if config_provider("sector", "solar_thermal")(w)
            else []
        ),
    output:
        snapshot_weightings=resources(
            "snapshot_weightings_base_s_{clusters}_elec_l{ll}_{opts}_{sector_opts}.csv"
        ),
    threads: 1
    resources:
        mem_mb=5000,
    log:
        logs("time_aggregation_base_s_{clusters}_elec_l{ll}_{opts}_{sector_opts}.log"),
    benchmark:
        benchmarks("time_aggregation_base_s_{clusters}_elec_l{ll}_{opts}_{sector_opts}")
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/time_aggregation.py"


def input_profile_offwind(w):
    return {
        f"profile_{tech}": resources("profile_{clusters}_" + tech + ".nc")
        for tech in ["offwind-ac", "offwind-dc", "offwind-float"]
        if (tech in config_provider("electricity", "renewable_carriers")(w))
    }


rule build_egs_potentials:
    params:
        snapshots=config_provider("snapshots"),
        sector=config_provider("sector"),
        costs=config_provider("costs"),
    input:
        egs_cost="data/egs_costs.json",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        air_temperature=(
            resources("temp_air_total_base_s_{clusters}.nc")
            if config_provider("sector", "enhanced_geothermal", "var_cf")
            else []
        ),
    output:
        egs_potentials=resources("egs_potentials_{clusters}.csv"),
        egs_overlap=resources("egs_overlap_{clusters}.csv"),
        egs_capacity_factors=resources("egs_capacity_factors_{clusters}.csv"),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs("build_egs_potentials_{clusters}.log"),
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/build_egs_potentials.py"


rule prepare_sector_network:
    params:
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
        co2_budget=config_provider("co2_budget"),
        conventional_carriers=config_provider(
            "existing_capacities", "conventional_carriers"
        ),
        foresight=config_provider("foresight"),
        costs=config_provider("costs"),
        sector=config_provider("sector"),
        industry=config_provider("industry"),
        renewable=config_provider("renewable"),
        lines=config_provider("lines"),
        pypsa_eur=config_provider("pypsa_eur"),
        length_factor=config_provider("lines", "length_factor"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        countries=config_provider("countries"),
        adjustments=config_provider("adjustments", "sector"),
        emissions_scope=config_provider("energy", "emissions"),
        biomass=config_provider("biomass"),
        RDIR=RDIR,
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        heat_systems=config_provider("sector", "heat_systems"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        unpack(input_profile_offwind),
        **rules.cluster_gas_network.output,
        **rules.build_gas_input_locations.output,
        snapshot_weightings=resources(
            "snapshot_weightings_base_s_{clusters}_elec_l{ll}_{opts}_{sector_opts}.csv"
        ),
        retro_cost=lambda w: (
            resources("retro_cost_base_s_{clusters}.csv")
            if config_provider("sector", "retrofitting", "retro_endogen")(w)
            else []
        ),
        floor_area=lambda w: (
            resources("floor_area_base_s_{clusters}.csv")
            if config_provider("sector", "retrofitting", "retro_endogen")(w)
            else []
        ),
        biomass_transport_costs=lambda w: (
            resources("biomass_transport_costs.csv")
            if config_provider("sector", "biomass_transport")(w)
            or config_provider("sector", "biomass_spatial")(w)
            else []
        ),
        sequestration_potential=lambda w: (
            resources("co2_sequestration_potential_base_s_{clusters}.csv")
            if config_provider(
                "sector", "regional_co2_sequestration_potential", "enable"
            )(w)
            else []
        ),
        network=resources("networks/base_s_{clusters}_elec_l{ll}_{opts}.nc"),
        eurostat="data/eurostat/Balances-April2023",
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        pop_weighted_heat_totals=resources("pop_weighted_heat_totals_s_{clusters}.csv"),
        shipping_demand=resources("shipping_demand_s_{clusters}.csv"),
        transport_demand=resources("transport_demand_s_{clusters}.csv"),
        transport_data=resources("transport_data_s_{clusters}.csv"),
        avail_profile=resources("avail_profile_s_{clusters}.csv"),
        dsm_profile=resources("dsm_profile_s_{clusters}.csv"),
        co2_totals_name=resources("co2_totals.csv"),
        co2="data/bundle/eea/UNFCCC_v23.csv",
        biomass_potentials=resources(
            "biomass_potentials_s_{clusters}_{planning_horizons}.csv"
        ),
        costs=lambda w: (
            resources("costs_{}.csv".format(config_provider("costs", "year")(w)))
            if config_provider("foresight")(w) == "overnight"
            else resources("costs_{planning_horizons}.csv")
        ),
        h2_cavern=resources("salt_cavern_potentials_s_{clusters}.csv"),
        busmap_s=resources("busmap_base_s.csv"),
        busmap=resources("busmap_base_s_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        industrial_demand=resources(
            "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
        ),
        hourly_heat_demand_total=resources(
            "hourly_heat_demand_total_base_s_{clusters}.nc"
        ),
        industrial_production=resources(
            "industrial_production_base_s_{clusters}_{planning_horizons}.csv"
        ),
        district_heat_share=resources(
            "district_heat_share_base_s_{clusters}_{planning_horizons}.csv"
        ),
        heating_efficiencies=resources("heating_efficiencies.csv"),
        temp_soil_total=resources("temp_soil_total_base_s_{clusters}.nc"),
        temp_air_total=resources("temp_air_total_base_s_{clusters}.nc"),
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
        solar_thermal_total=lambda w: (
            resources("solar_thermal_total_base_s_{clusters}.nc")
            if config_provider("sector", "solar_thermal")(w)
            else []
        ),
        egs_potentials=lambda w: (
            resources("egs_potentials_{clusters}.csv")
            if config_provider("sector", "enhanced_geothermal", "enable")(w)
            else []
        ),
        egs_overlap=lambda w: (
            resources("egs_overlap_{clusters}.csv")
            if config_provider("sector", "enhanced_geothermal", "enable")(w)
            else []
        ),
        egs_capacity_factors=lambda w: (
            resources("egs_capacity_factors_{clusters}.csv")
            if config_provider("sector", "enhanced_geothermal", "enable")(w)
            else []
        ),
    output:
        RESULTS
        + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        RESULTS
        + "logs/prepare_sector_network_base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/prepare_sector_network/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/prepare_sector_network.py"
