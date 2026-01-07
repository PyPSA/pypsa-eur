# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT


rule build_population_layouts:
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        urban_percent=rules.retrieve_worldbank_urban_population.output["csv"],
        cutout=lambda w: input_cutout(w),
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
    script:
        "../scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=lambda w: input_cutout(w),
    output:
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
    log:
        logs("build_clustered_population_layouts_s_{clusters}.log"),
    resources:
        mem_mb=10000,
    benchmark:
        benchmarks("build_clustered_population_layouts/s_{clusters}")
    script:
        "../scripts/build_clustered_population_layouts.py"


rule build_clustered_solar_rooftop_potentials:
    input:
        pop_layout=resources("pop_layout_total.nc"),
        class_regions=resources("regions_by_class_{clusters}_solar.geojson"),
        cutout=lambda w: input_cutout(w),
    output:
        potentials=resources("solar_rooftop_potentials_s_{clusters}.csv"),
    log:
        logs("build_clustered_solar_rooftop_potentials_s_{clusters}.log"),
    resources:
        mem_mb=10000,
    benchmark:
        benchmarks("build_clustered_solar_rooftop_potentials/s_{clusters}")
    script:
        "../scripts/build_clustered_solar_rooftop_potentials.py"


rule build_simplified_population_layouts:
    input:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
        regions_onshore=resources("regions_onshore_base_s.geojson"),
        cutout=lambda w: input_cutout(w),
    output:
        clustered_pop_layout=resources("pop_layout_base_s.csv"),
    resources:
        mem_mb=10000,
    log:
        logs("build_simplified_population_layouts_s"),
    benchmark:
        benchmarks("build_simplified_population_layouts/s")
    script:
        "../scripts/build_clustered_population_layouts.py"


rule build_gas_network:
    input:
        gas_network=rules.retrieve_gas_infrastructure_data.output["gas_network"],
    output:
        cleaned_gas_network=resources("gas_network.csv"),
    resources:
        mem_mb=4000,
    log:
        logs("build_gas_network.log"),
    benchmark:
        benchmarks("build_gas_network")
    script:
        "../scripts/build_gas_network.py"


rule build_gas_input_locations:
    input:
        gem="data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        entry=rules.retrieve_gas_infrastructure_data.output["entry"],
        storage=rules.retrieve_gas_infrastructure_data.output["storage"],
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
    benchmark:
        benchmarks("build_gas_input_locations/s_{clusters}")
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
    benchmark:
        benchmarks("cluster_gas_network/s_{clusters}")
    script:
        "../scripts/cluster_gas_network.py"


rule build_daily_heat_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        pop_layout=resources("pop_layout_total.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=lambda w: input_cutout(
            w, config_provider("sector", "heat_demand_cutout")(w)
        ),
    output:
        heat_demand=resources("daily_heat_demand_total_base_s_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 8
    log:
        logs("build_daily_heat_demand_total_s_{clusters}.loc"),
    benchmark:
        benchmarks("build_daily_heat_demand/total_s_{clusters}")
    script:
        "../scripts/build_daily_heat_demand.py"


rule build_hourly_heat_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        sector=config_provider("sector"),
    input:
        heat_profile="data/heat_load_profile_BDEW.csv",
        heat_demand=resources("daily_heat_demand_total_base_s_{clusters}.nc"),
    output:
        heat_demand=resources("hourly_heat_demand_total_base_s_{clusters}.nc"),
        heat_dsm_profile=resources(
            "residential_heat_dsm_profile_total_base_s_{clusters}.csv"
        ),
    resources:
        mem_mb=2000,
    threads: 8
    log:
        logs("build_hourly_heat_demand_total_s_{clusters}.loc"),
    benchmark:
        benchmarks("build_hourly_heat_demand/total_s_{clusters}")
    script:
        "../scripts/build_hourly_heat_demand.py"


rule build_temperature_profiles:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        pop_layout=resources("pop_layout_total.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=lambda w: input_cutout(
            w, config_provider("sector", "heat_demand_cutout")(w)
        ),
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
        drop_leap_day=config_provider("enable", "drop_leap_day"),
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
    script:
        "../scripts/build_central_heating_temperature_profiles/run.py"


rule build_dh_areas:
    params:
        handle_missing_countries=config_provider(
            "sector", "district_heating", "dh_areas", "handle_missing_countries"
        ),
        countries=config_provider("countries"),
    input:
        dh_areas=rules.retrieve_dh_areas.output["dh_areas"],
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        dh_areas=resources("dh_areas_base_s_{clusters}.geojson"),
    resources:
        mem_mb=2000,
    log:
        logs("build_dh_areas_s_{clusters}.log"),
    benchmark:
        benchmarks("build_dh_areas_s/s_{clusters}")
    script:
        "../scripts/build_dh_areas.py"


rule build_geothermal_heat_potential:
    params:
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        countries=config_provider("countries"),
        constant_temperature_celsius=config_provider(
            "sector",
            "district_heating",
            "limited_heat_sources",
            "geothermal",
            "constant_temperature_celsius",
        ),
        ignore_missing_regions=config_provider(
            "sector",
            "district_heating",
            "limited_heat_sources",
            "geothermal",
            "ignore_missing_regions",
        ),
    input:
        isi_heat_potentials=rules.retrieve_geothermal_heat_utilisation_potentials.output[
            "isi_heat_potentials"
        ],
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        lau_regions=rules.retrieve_lau_regions.output["zip"],
    output:
        heat_source_power=resources(
            "heat_source_power_geothermal_base_s_{clusters}.csv"
        ),
    resources:
        mem_mb=2000,
    log:
        logs("build_heat_source_potentials_geothermal_s_{clusters}.log"),
    benchmark:
        benchmarks("build_heat_source_potentials/geothermal_s_{clusters}")
    script:
        "../scripts/build_geothermal_heat_potential.py"


rule build_ates_potentials:
    params:
        max_top_temperature=config_provider(
            "sector",
            "district_heating",
            "ates",
            "max_top_temperature",
        ),
        min_bottom_temperature=config_provider(
            "sector",
            "district_heating",
            "ates",
            "min_bottom_temperature",
        ),
        suitable_aquifer_types=config_provider(
            "sector",
            "district_heating",
            "ates",
            "suitable_aquifer_types",
        ),
        aquifer_volumetric_heat_capacity=config_provider(
            "sector",
            "district_heating",
            "ates",
            "aquifer_volumetric_heat_capacity",
        ),
        fraction_of_aquifer_area_available=config_provider(
            "sector",
            "district_heating",
            "ates",
            "fraction_of_aquifer_area_available",
        ),
        effective_screen_length=config_provider(
            "sector",
            "district_heating",
            "ates",
            "effective_screen_length",
        ),
        dh_area_buffer=config_provider(
            "sector",
            "district_heating",
            "dh_areas",
            "buffer",
        ),
        ignore_missing_regions=config_provider(
            "sector",
            "district_heating",
            "ates",
            "ignore_missing_regions",
        ),
        countries=config_provider("countries"),
    input:
        aquifer_shapes_shp=rules.retrieve_aquifer_data_bgr.output["aquifer_shapes"][0],
        dh_areas=resources("dh_areas_base_s_{clusters}.geojson"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
    output:
        ates_potentials=resources(
            "ates_potentials_base_s_{clusters}_{planning_horizons}.csv"
        ),
    resources:
        mem_mb=2000,
    log:
        logs("build_ates_potentials_s_{clusters}_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_ates_potentials_geothermal_s_{clusters}_{planning_horizons}")
    script:
        "../scripts/build_ates_potentials.py"


def input_hera_data(w) -> dict[str, str]:
    """
    Generate input file paths for HERA river discharge and ambient temperature data.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Snakemake wildcards object.

    Returns
    -------
    dict[str, str]
        Dictionary mapping keys like "hera_river_discharge_{year}" and
        "hera_ambient_temperature_{year}" to NetCDF file paths.
    """
    if config_provider("atlite", "default_cutout")(w) == "be-03-2013-era5":
        hera_data_key = "be_2013-03-01_to_2013-03-08"
        return {
            "hera_river_discharge_2013": f"data/hera_{hera_data_key}/river_discharge_{hera_data_key}.nc",
            "hera_ambient_temperature_2013": f"data/hera_{hera_data_key}/ambient_temp_{hera_data_key}.nc",
        }
    else:
        from scripts._helpers import get_snapshots

        # Get all snapshots and extract unique years
        snapshots_config = config_provider("snapshots")(w)
        snapshots = get_snapshots(snapshots_config)
        unique_years = snapshots.year.unique()

        # Create dictionary with year-specific keys
        result = {}
        for year in unique_years:
            result[f"hera_river_discharge_{year}"] = (
                f"data/hera_{year}/river_discharge_{year}.nc"
            )
            result[f"hera_ambient_temperature_{year}"] = (
                f"data/hera_{year}/ambient_temp_{year}.nc"
            )

        return result


rule build_river_heat_potential:
    params:
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        snapshots=config_provider("snapshots"),
        dh_area_buffer=config_provider(
            "sector", "district_heating", "dh_areas", "buffer"
        ),
        enable_heat_source_maps=config_provider("plotting", "enable_heat_source_maps"),
    input:
        unpack(input_hera_data),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        dh_areas=resources("dh_areas_base_s_{clusters}.geojson"),
    output:
        heat_source_power=resources(
            "heat_source_power_river_water_base_s_{clusters}.csv"
        ),
        heat_source_temperature=resources("temp_river_water_base_s_{clusters}.nc"),
        heat_source_temperature_temporal_aggregate=resources(
            "temp_river_water_base_s_{clusters}_temporal_aggregate.nc"
        ),
        heat_source_energy_temporal_aggregate=resources(
            "heat_source_energy_river_water_base_s_{clusters}_temporal_aggregate.nc"
        ),
    resources:
        mem_mb=20000,
    log:
        logs("build_river_water_heat_potential_base_s_{clusters}.log"),
    benchmark:
        benchmarks("build_river_water_heat_potential_base_s_{clusters}")
    threads: 1
    script:
        "../scripts/build_surface_water_heat_potentials/build_river_water_heat_potential.py"


def input_heat_source_temperature(
    w,
    replace_names: dict[str, str] = {
        "air": "air_total",
        "ground": "soil_total",
        "ptes": "ptes_top_profiles",
    },
) -> dict[str, str]:
    """
    Generate input file paths for heat source temperature profiles.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Snakemake wildcards object.
    replace_names : dict[str, str], optional
        Mapping to transform heat source names to file naming conventions.

    Returns
    -------
    dict[str, str]
        Dictionary mapping keys like "temp_{heat_source_name}" to NetCDF file paths
        for heat sources that require temperature profiles (excludes constant
        temperature sources).
    """

    heat_pump_sources = set(
        config_provider("sector", "heat_pump_sources", "urban central")(w)
    ).union(
        config_provider("sector", "heat_pump_sources", "urban decentral")(w),
        config_provider("sector", "heat_pump_sources", "rural")(w),
    )

    is_limited_heat_source = {
        heat_source_name: heat_source_name
        in config_provider("sector", "district_heating", "limited_heat_sources")(w)
        for heat_source_name in heat_pump_sources
    }

    has_constant_temperature = {
        heat_source_name: (
            False
            if not is_limited_heat_source[heat_source_name]
            else config_provider(
                "sector",
                "district_heating",
                "limited_heat_sources",
                heat_source_name,
                "constant_temperature_celsius",
            )(w)
        )
        for heat_source_name in heat_pump_sources
    }

    # replace names for soil and air temperature files
    return {
        f"temp_{heat_source_name}": resources(
            "temp_"
            + replace_names.get(heat_source_name, heat_source_name)
            + "_base_s_{clusters}"
            + ("_{planning_horizons}" if heat_source_name == "ptes" else "")
            + ".nc"
        )
        for heat_source_name in heat_pump_sources
        # remove heat sources with constant temperature - i.e. no temperature profile file
        if not has_constant_temperature[heat_source_name]
    }


def input_seawater_temperature(w) -> dict[str, str]:
    """
    Generate input file paths for seawater temperature data.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Snakemake wildcards object.

    Returns
    -------
    dict[str, str]
        Dictionary mapping keys like "seawater_temperature_{year}" to NetCDF file paths.
    """

    # Import here to avoid circular imports
    from scripts._helpers import get_snapshots

    # Get all snapshots and extract unique years
    snapshots_config = config_provider("snapshots")(w)
    snapshots = get_snapshots(snapshots_config)
    unique_years = snapshots.year.unique()

    # Create dictionary with year-specific keys
    return {
        f"seawater_temperature_{year}": f"data/seawater_temperature_{year}.nc"
        for year in unique_years
    }


rule build_sea_heat_potential:
    params:
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        snapshots=config_provider("snapshots"),
        dh_area_buffer=config_provider(
            "sector", "district_heating", "dh_areas", "buffer"
        ),
    input:
        # seawater_temperature=lambda w: input_seawater_temperature(w),
        unpack(input_seawater_temperature),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        dh_areas=resources("dh_areas_base_s_{clusters}.geojson"),
    output:
        heat_source_temperature=resources("temp_sea_water_base_s_{clusters}.nc"),
        heat_source_temperature_temporal_aggregate=resources(
            "temp_sea_water_base_s_{clusters}_temporal_aggregate.nc"
        ),
    resources:
        mem_mb=10000,
    log:
        logs("build_sea_water_heat_potential_base_s_{clusters}.log"),
    benchmark:
        benchmarks("build_sea_water_heat_potential_base_s_{clusters}")
    threads: config["atlite"].get("nprocesses", 4)
    script:
        "../scripts/build_surface_water_heat_potentials/build_sea_water_heat_potential.py"


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
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        snapshots=config_provider("snapshots"),
    input:
        unpack(input_heat_source_temperature),
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        cop_profiles=resources("cop_profiles_base_s_{clusters}_{planning_horizons}.nc"),
    resources:
        mem_mb=20000,
    log:
        logs("build_cop_profiles_s_{clusters}_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_cop_profiles/s_{clusters}_{planning_horizons}")
    script:
        "../scripts/build_cop_profiles/run.py"


rule build_ptes_operations:
    params:
        max_ptes_top_temperature=config_provider(
            "sector",
            "district_heating",
            "ptes",
            "max_top_temperature",
        ),
        min_ptes_bottom_temperature=config_provider(
            "sector",
            "district_heating",
            "ptes",
            "min_bottom_temperature",
        ),
        snapshots=config_provider("snapshots"),
    input:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
    output:
        ptes_direct_utilisation_profiles=resources(
            "ptes_direct_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        ptes_top_temperature_profiles=resources(
            "temp_ptes_top_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        ptes_e_max_pu_profiles=resources(
            "ptes_e_max_pu_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
    resources:
        mem_mb=2000,
    log:
        logs("build_ptes_operations_s_{clusters}_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_ptes_operations_s_{clusters}_{planning_horizons}")
    script:
        "../scripts/build_ptes_operations/run.py"


rule build_direct_heat_source_utilisation_profiles:
    params:
        direct_utilisation_heat_sources=config_provider(
            "sector", "district_heating", "direct_utilisation_heat_sources"
        ),
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        snapshots=config_provider("snapshots"),
    input:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
    output:
        direct_heat_source_utilisation_profiles=resources(
            "direct_heat_source_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
    resources:
        mem_mb=20000,
    log:
        logs(
            "build_direct_heat_source_utilisation_profiles_s_{clusters}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "build_direct_heat_source_utilisation_profiles/s_{clusters}_{planning_horizons}"
        )
    script:
        "../scripts/build_direct_heat_source_utilisation_profiles.py"


rule build_solar_thermal_profiles:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        solar_thermal=config_provider("solar_thermal"),
    input:
        pop_layout=resources("pop_layout_total.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        cutout=lambda w: input_cutout(w, config_provider("solar_thermal", "cutout")(w)),
    output:
        solar_thermal=resources("solar_thermal_total_base_s_{clusters}.nc"),
    resources:
        mem_mb=20000,
    threads: 16
    log:
        logs("build_solar_thermal_profiles_total_s_{clusters}.log"),
    benchmark:
        benchmarks("build_solar_thermal_profiles/total_{clusters}")
    script:
        "../scripts/build_solar_thermal_profiles.py"


rule build_energy_totals:
    params:
        countries=config_provider("countries"),
        energy=config_provider("energy"),
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        co2=rules.retrieve_ghg_emissions.output["csv"],
        swiss="data/switzerland-new_format-all_years.csv",
        swiss_transport=f"{BFS_ROAD_VEHICLE_STOCK_DATASET['folder']}/vehicle_stock.csv",
        idees=rules.retrieve_jrc_idees.output["directory"],
        district_heat_share="data/district_heat_share.csv",
        eurostat=rules.retrieve_eurostat_balances.output["directory"],
        eurostat_households=rules.retrieve_eurostat_household_balances.output["csv"],
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
    script:
        "../scripts/build_energy_totals.py"


if (COUNTRY_HDD_DATASET := dataset_version("country_hdd"))["source"] in ["build"]:

    # This rule uses one or multiple cutouts.
    # To update the output files to include a new year, e.g. 2025 using an existing cutout,
    # either create a new cutout covering the whole timespan or add another cutout that covers the additional year(s).
    # E.g. cutouts=[<cutout for 1940-2024>, <cutout for 2025-2025>]
    rule build_country_hdd:
        input:
            cutouts=["cutouts/europe-1940-2024-era5.nc"],
            country_shapes=resources("country_shapes.geojson"),
        output:
            era5_hdd=f"{COUNTRY_HDD_DATASET["folder"]}/era5-HDD-per-country.csv",
        log:
            logs("build_country_hdd.log"),
        benchmark:
            benchmarks("build_country_hdd")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_country_hdd.py"


rule build_heat_totals:
    input:
        hdd=f"{COUNTRY_HDD_DATASET["folder"]}/era5-HDD-per-country.csv",
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
    script:
        "../scripts/build_heat_totals.py"


rule build_biomass_potentials:
    params:
        biomass=config_provider("biomass"),
    input:
        enspreso_biomass=rules.retrieve_enspreso_biomass.output["xlsx"],
        eurostat=rules.retrieve_eurostat_balances.output["directory"],
        nuts2=rules.retrieve_eu_nuts_2013.output["shapes_level_2"],
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        nuts3_population=ancient(rules.retrieve_nuts3_population.output["gz"]),
        swiss_cantons=ancient("data/ch_cantons.csv"),
        swiss_population=rules.retrieve_bfs_gdp_and_population.output["xlsx"],
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
        mem_mb=2000,
    log:
        logs("build_biomass_potentials_s_{clusters}_{planning_horizons}.log"),
    benchmark:
        benchmarks("build_biomass_potentials_s_{clusters}_{planning_horizons}")
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
    script:
        "../scripts/build_biomass_transport_costs.py"


rule build_co2_sequestration_potentials:
    input:
        storage_table=rules.retrieve_co2stop.output["storage_table"],
        storage_map=rules.retrieve_co2stop.output["storage_map"],
        traps_table1=rules.retrieve_co2stop.output["traps_table1"],
        traps_table2=rules.retrieve_co2stop.output["traps_table2"],
        traps_table3=rules.retrieve_co2stop.output["traps_table3"],
        traps_map=rules.retrieve_co2stop.output["traps_map"],
    output:
        resources("co2_sequestration_potentials.geojson"),
    threads: 1
    resources:
        mem_mb=4000,
    log:
        logs("build_co2_sequestration_potentials.log"),
    benchmark:
        benchmarks("build_co2_sequestration_potentials")
    script:
        "../scripts/build_co2_sequestration_potentials.py"


rule build_clustered_co2_sequestration_potentials:
    params:
        sequestration_potential=config_provider(
            "sector", "regional_co2_sequestration_potential"
        ),
    input:
        sequestration_potential=resources("co2_sequestration_potentials.geojson"),
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
        logs("build_clustered_co2_sequestration_potentials_{clusters}.log"),
    benchmark:
        benchmarks("build_clustered_co2_sequestration_potentials_{clusters}")
    script:
        "../scripts/build_clustered_co2_sequestration_potentials.py"


rule build_salt_cavern_potentials:
    input:
        salt_caverns=rules.retrieve_h2_salt_caverns.output["geojson"],
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
    script:
        "../scripts/build_salt_cavern_potentials.py"


rule build_ammonia_production:
    input:
        usgs=rules.retrieve_nitrogen_statistics.output["xlsx"],
    output:
        ammonia_production=resources("ammonia_production.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_ammonia_production.log"),
    benchmark:
        benchmarks("build_ammonia_production")
    script:
        "../scripts/build_ammonia_production.py"


rule build_industry_sector_ratios:
    params:
        industry=config_provider("industry"),
        ammonia=config_provider("sector", "ammonia", default=False),
    input:
        ammonia_production=resources("ammonia_production.csv"),
        idees=rules.retrieve_jrc_idees.output["directory"],
    output:
        industry_sector_ratios=resources("industry_sector_ratios.csv"),
    threads: 1
    resources:
        mem_mb=1000,
    log:
        logs("build_industry_sector_ratios.log"),
    benchmark:
        benchmarks("build_industry_sector_ratios")
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
    script:
        "../scripts/build_industry_sector_ratios_intermediate.py"


rule build_industrial_production_per_country:
    params:
        industry=config_provider("industry"),
        countries=config_provider("countries"),
    input:
        ch_industrial_production="data/ch_industrial_production_per_subsector.csv",
        ammonia_production=resources("ammonia_production.csv"),
        eurostat=rules.retrieve_eurostat_balances.output["directory"],
        jrc=rules.retrieve_jrc_idees.output["directory"],
    output:
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    threads: 8
    resources:
        mem_mb=2000,
    log:
        logs("build_industrial_production_per_country.log"),
    benchmark:
        benchmarks("build_industrial_production_per_country")
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
        hotmaps=rules.retrieve_hotmaps_industrial_sites.output["csv"],
        gem_gspt=rules.retrieve_gem_steel_plant_tracker.output["xlsx"],
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
    script:
        "../scripts/build_industrial_energy_demand_per_node.py"


rule build_industrial_energy_demand_per_country_today:
    params:
        countries=config_provider("countries"),
        industry=config_provider("industry"),
        ammonia=config_provider("sector", "ammonia", default=False),
    input:
        transformation_output_coke=resources("transformation_output_coke.csv"),
        jrc=rules.retrieve_jrc_idees.output["directory"],
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    output:
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
    threads: 8
    resources:
        mem_mb=2000,
    log:
        logs("build_industrial_energy_demand_per_country_today.log"),
    benchmark:
        benchmarks("build_industrial_energy_demand_per_country_today")
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
    script:
        "../scripts/build_retro_cost.py"


rule build_population_weighted_energy_totals:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
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
    benchmark:
        benchmarks("build_population_weighted_{kind}_totals_{clusters}")
    script:
        "../scripts/build_population_weighted_energy_totals.py"


rule build_shipping_demand:
    input:
        ports=rules.retrieve_attributed_ports.output["json"],
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
    benchmark:
        benchmarks("build_shipping_demand/s_{clusters}")
    script:
        "../scripts/build_shipping_demand.py"


if MOBILITY_PROFILES_DATASET["source"] in ["build"]:

    rule build_mobility_profiles:
        params:
            sector=config_provider("sector"),
        input:
            zip_files=storage(
                expand(
                    MOBILITY_PROFILES_DATASET["url"],
                    year=[2010, 2011, 2012, 2013, 2014],
                    street_type=["A", "B"],
                ),
            ),
        output:
            raw_files=directory(MOBILITY_PROFILES_DATASET["folder"] / "raw"),
            kfz=MOBILITY_PROFILES_DATASET["folder"] / "kfz.csv",
            pkw=MOBILITY_PROFILES_DATASET["folder"] / "pkw.csv",
        threads: 1
        resources:
            mem_mb=5000,
        log:
            logs("build_mobility_profiles.log"),
        benchmark:
            benchmarks("build_mobility_profiles")
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/build_mobility_profiles.py"


rule build_transport_demand:
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        sector=config_provider("sector"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    input:
        network=resources("networks/base_s.nc"),
        clustered_pop_layout=resources("pop_layout_base_s_{clusters}.csv"),
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        transport_data=resources("transport_data.csv"),
        traffic_data_KFZ=f"{MOBILITY_PROFILES_DATASET["folder"]}/kfz.csv",
        traffic_data_Pkw=f"{MOBILITY_PROFILES_DATASET["folder"]}/pkw.csv",
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
    benchmark:
        benchmarks("build_transport_demand/s_{clusters}")
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
    benchmark:
        benchmarks("build_district_heat_share_{clusters}_{planning_horizons}")
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
    script:
        "../scripts/build_existing_heating_distribution.py"


rule time_aggregation:
    params:
        time_resolution=config_provider("clustering", "temporal"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        solver_name=config_provider("solving", "solver", "name"),
    input:
        network=resources("networks/base_s_{clusters}_elec_{opts}.nc"),
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
            "snapshot_weightings_base_s_{clusters}_elec_{opts}_{sector_opts}.csv"
        ),
    threads: 1
    resources:
        mem_mb=5000,
    log:
        logs("time_aggregation_base_s_{clusters}_elec_{opts}_{sector_opts}.log"),
    benchmark:
        benchmarks("time_aggregation_base_s_{clusters}_elec_{opts}_{sector_opts}")
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
        drop_leap_day=config_provider("enable", "drop_leap_day"),
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
    benchmark:
        benchmarks("build_egs_potentials_{clusters}")
    script:
        "../scripts/build_egs_potentials.py"


def input_heat_source_power(w):

    return {
        heat_source_name: resources(
            "heat_source_power_" + heat_source_name + "_base_s_{clusters}.csv"
        )
        for heat_source_name in config_provider(
            "sector", "heat_pump_sources", "urban central"
        )(w)
        if heat_source_name
        in config_provider("sector", "district_heating", "limited_heat_sources")(
            w
        ).keys()
    }


rule prepare_sector_network:
    params:
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
        co2_budget=config_provider("co2_budget"),
        conventional_carriers=config_provider(
            "existing_capacities", "conventional_carriers"
        ),
        foresight=config_provider("foresight"),
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
        emission_prices=config_provider("costs", "emission_prices"),
        electricity=config_provider("electricity"),
        biomass=config_provider("biomass"),
        RDIR=RDIR,
        heat_pump_sources=config_provider("sector", "heat_pump_sources"),
        heat_systems=config_provider("sector", "heat_systems"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
        direct_utilisation_heat_sources=config_provider(
            "sector", "district_heating", "direct_utilisation_heat_sources"
        ),
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        temperature_limited_stores=config_provider(
            "sector", "district_heating", "temperature_limited_stores"
        ),
    input:
        unpack(input_profile_offwind),
        unpack(input_heat_source_power),
        **rules.cluster_gas_network.output,
        **rules.build_gas_input_locations.output,
        snapshot_weightings=resources(
            "snapshot_weightings_base_s_{clusters}_elec_{opts}_{sector_opts}.csv"
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
        network=resources("networks/base_s_{clusters}_elec_{opts}.nc"),
        eurostat=rules.retrieve_eurostat_balances.output["directory"],
        pop_weighted_energy_totals=resources(
            "pop_weighted_energy_totals_s_{clusters}.csv"
        ),
        pop_weighted_heat_totals=resources("pop_weighted_heat_totals_s_{clusters}.csv"),
        shipping_demand=resources("shipping_demand_s_{clusters}.csv"),
        transport_demand=resources("transport_demand_s_{clusters}.csv"),
        transport_data=resources("transport_data_s_{clusters}.csv"),
        avail_profile=resources("avail_profile_s_{clusters}.csv"),
        dsm_profile=resources("dsm_profile_s_{clusters}.csv"),
        heat_dsm_profile=resources(
            "residential_heat_dsm_profile_total_base_s_{clusters}.csv"
        ),
        co2_totals_name=resources("co2_totals.csv"),
        co2=rules.retrieve_ghg_emissions.output["csv"],
        biomass_potentials=resources(
            "biomass_potentials_s_{clusters}_{planning_horizons}.csv"
        ),
        costs=lambda w: (
            resources(f"costs_{config_provider("costs", "year")(w)}_processed.csv")
            if config_provider("foresight")(w) == "overnight"
            else resources("costs_{planning_horizons}_processed.csv")
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
        ptes_e_max_pu_profiles=lambda w: (
            resources(
                "ptes_e_max_pu_profiles_base_s_{clusters}_{planning_horizons}.nc"
            )
            if config_provider(
                "sector", "district_heating", "ptes", "dynamic_capacity"
            )(w)
            else []
        ),
        ptes_direct_utilisation_profiles=lambda w: (
            resources(
                "ptes_direct_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc"
            )
            if config_provider(
                "sector", "district_heating", "ptes", "supplemental_heating", "enable"
            )(w)
            else []
        ),
        solar_thermal_total=lambda w: (
            resources("solar_thermal_total_base_s_{clusters}.nc")
            if config_provider("sector", "solar_thermal")(w)
            else []
        ),
        solar_rooftop_potentials=lambda w: (
            resources("solar_rooftop_potentials_s_{clusters}.csv")
            if "solar" in config_provider("electricity", "renewable_carriers")(w)
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
        direct_heat_source_utilisation_profiles=resources(
            "direct_heat_source_utilisation_profiles_base_s_{clusters}_{planning_horizons}.nc"
        ),
        ates_potentials=lambda w: (
            resources("ates_potentials_base_s_{clusters}_{planning_horizons}.csv")
            if config_provider("sector", "district_heating", "ates", "enable")(w)
            else []
        ),
    output:
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
    threads: 1
    resources:
        mem_mb=2000,
    log:
        logs(
            "prepare_sector_network_base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.log"
        ),
    benchmark:
        benchmarks(
            "prepare_sector_network/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}"
        )
    script:
        "../scripts/prepare_sector_network.py"
