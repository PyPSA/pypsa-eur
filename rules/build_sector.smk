# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule build_population_layouts:
    """Maps population (total, urban, rural) onto weather cutout grid cells from NUTS3 shapes."""
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
    benchmark:
        benchmarks("build_population_layouts")
    threads: 8
    resources:
        mem_mb=20000,
    script:
        scripts("build_population_layouts.py")


rule build_clustered_population_layouts:
    """Aggregates total, urban and rural population layouts to clustered model regions."""
    input:
        pop_layout_total=resources("pop_layout_total.nc"),
        pop_layout_urban=resources("pop_layout_urban.nc"),
        pop_layout_rural=resources("pop_layout_rural.nc"),
        onshore_regions=resources("onshore_regions.geojson"),
        cutout=lambda w: input_cutout(w),
    output:
        clustered_pop_layout=resources("pop_layout.csv"),
    log:
        logs("build_clustered_population_layouts.log"),
    benchmark:
        benchmarks("build_clustered_population_layouts")
    resources:
        mem_mb=10000,
    script:
        scripts("build_clustered_population_layouts.py")


rule build_solar_rooftop_potentials:
    """Computes solar rooftop potentials per resource class for all clustered model regions."""
    input:
        pop_layout=resources("pop_layout_total.nc"),
        class_regions=resources("regions_by_class_solar.geojson"),
        cutout=lambda w: input_cutout(w),
    output:
        potentials=resources("solar_rooftop_potentials.csv"),
    log:
        logs("build_solar_rooftop_potentials.log"),
    benchmark:
        benchmarks("build_solar_rooftop_potentials")
    resources:
        mem_mb=10000,
    script:
        scripts("build_solar_rooftop_potentials.py")


rule build_gas_network:
    """Preprocesses the SciGRID_gas transmission network into cleaned pipeline segments."""
    input:
        gas_network=rules.retrieve_gas_infrastructure_data.output["gas_network"],
    output:
        cleaned_gas_network=resources("gas_network.csv"),
    log:
        logs("build_gas_network.log"),
    benchmark:
        benchmarks("build_gas_network")
    resources:
        mem_mb=4000,
    script:
        scripts("build_gas_network.py")


rule build_gas_input_locations:
    """Builds fossil gas import locations from entry points, LNG terminals and production sites."""
    input:
        gem="data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        entry=rules.retrieve_gas_infrastructure_data.output["entry"],
        storage=rules.retrieve_gas_infrastructure_data.output["storage"],
        onshore_regions=resources("onshore_regions.geojson"),
        offshore_regions=resources("offshore_regions.geojson"),
    output:
        gas_input_nodes=resources("gas_input_locations.geojson"),
        gas_input_nodes_simplified=resources("gas_input_locations_simplified.csv"),
    log:
        logs("build_gas_input_locations.log"),
    benchmark:
        benchmarks("build_gas_input_locations")
    resources:
        mem_mb=2000,
    script:
        scripts("build_gas_input_locations.py")


rule cluster_gas_network:
    """Clusters the gas transmission network pipelines to clustered model regions."""
    input:
        cleaned_gas_network=resources("gas_network.csv"),
        onshore_regions=resources("onshore_regions.geojson"),
        offshore_regions=resources("offshore_regions.geojson"),
    output:
        clustered_gas_network=resources("gas_network_clustered.csv"),
    log:
        logs("cluster_gas_network.log"),
    benchmark:
        benchmarks("cluster_gas_network")
    resources:
        mem_mb=4000,
    script:
        scripts("cluster_gas_network.py")


rule build_daily_heat_demand:
    """Builds daily heat demand time series per region from cutout temperatures and population."""
    input:
        pop_layout=resources("pop_layout_total.nc"),
        onshore_regions=resources("onshore_regions.geojson"),
        cutout=lambda w: input_cutout(
            w, config_provider("sector", "heat_demand_cutout")(w)
        ),
    output:
        heat_demand=resources("daily_heat_demand_total.nc"),
    log:
        logs("build_daily_heat_demand_total.log"),
    benchmark:
        benchmarks("build_daily_heat_demand_total")
    threads: 8
    resources:
        mem_mb=20000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    script:
        scripts("build_daily_heat_demand.py")


rule build_hourly_heat_demand:
    """Disaggregates daily heat demand into hourly profiles with standard load curves."""
    input:
        heat_profile="data/heat_load_profile_BDEW.csv",
        heat_demand=resources("daily_heat_demand_total.nc"),
    output:
        heat_demand=resources("hourly_heat_demand_total.nc"),
        heat_dsm_profile=resources("residential_heat_dsm_profile.csv"),
    log:
        logs("build_hourly_heat_demand_total.loc"),
    benchmark:
        benchmarks("build_hourly_heat_demand_total")
    threads: 8
    resources:
        mem_mb=2000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        sector=config_provider("sector"),
    script:
        scripts("build_hourly_heat_demand.py")


rule build_temperature_profiles:
    """Builds population-weighted air and soil temperature time series per clustered region."""
    input:
        pop_layout=resources("pop_layout_total.nc"),
        onshore_regions=resources("onshore_regions.geojson"),
        cutout=lambda w: input_cutout(
            w, config_provider("sector", "heat_demand_cutout")(w)
        ),
    output:
        temp_soil=resources("temp_soil_total.nc"),
        temp_air=resources("temp_air_total.nc"),
    log:
        logs("build_temperature_profiles_total.log"),
    benchmark:
        benchmarks("build_temperature_profiles/total")
    threads: 8
    resources:
        mem_mb=20000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    script:
        scripts("build_temperature_profiles.py")


rule build_central_heating_temperature_profiles:
    """Approximates district heating forward and return temperature profiles from air temperature."""
    input:
        temp_air_total=resources("temp_air_total.nc"),
        onshore_regions=resources("onshore_regions.geojson"),
    output:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_{horizon}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_{horizon}.nc"
        ),
    log:
        logs("build_central_heating_temperature_profiles_{horizon}.log"),
    benchmark:
        benchmarks("build_central_heating_temperature_profiles_{horizon}")
    resources:
        mem_mb=20000,
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
    script:
        scripts("build_central_heating_temperature_profiles/run.py")


rule build_dh_areas:
    """Builds district heating area shapes and fills in countries missing from the source data."""
    input:
        dh_areas=rules.retrieve_dh_areas.output["dh_areas"],
        onshore_regions=resources("onshore_regions.geojson"),
    output:
        dh_areas=resources("dh_areas.geojson"),
    log:
        logs("build_dh_areas.log"),
    benchmark:
        benchmarks("build_dh_areas")
    resources:
        mem_mb=2000,
    script:
        scripts("build_dh_areas.py")


rule build_geothermal_heat_potential:
    """Aggregates LAU-level geothermal heat potentials to technical potentials per model region."""
    input:
        isi_heat_potentials=rules.retrieve_geothermal_heat_utilisation_potentials.output[
            "isi_heat_potentials"
        ],
        onshore_regions=resources("onshore_regions.geojson"),
        lau_regions=rules.retrieve_lau_regions.output["zip"],
    output:
        heat_source_power=resources("heat_source_power_geothermal.csv"),
    log:
        logs("build_heat_source_potentials_geothermal.log"),
    benchmark:
        benchmarks("build_heat_source_potentials/geothermal")
    resources:
        mem_mb=2000,
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
    script:
        scripts("build_geothermal_heat_potential.py")


rule build_ates_potentials:
    """Computes aquifer thermal energy storage potentials per region from aquifers and heating areas."""
    input:
        aquifer_shapes_shp=rules.retrieve_aquifer_data_bgr.output["aquifer_shapes"][0],
        dh_areas=resources("dh_areas.geojson"),
        onshore_regions=resources("onshore_regions.geojson"),
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_{horizon}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_{horizon}.nc"
        ),
    output:
        ates_potentials=resources("ates_potentials_{horizon}.csv"),
    log:
        logs("build_ates_potentials_{horizon}.log"),
    benchmark:
        benchmarks("build_ates_potentials_geothermal_{horizon}")
    resources:
        mem_mb=2000,
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
    script:
        scripts("build_ates_potentials.py")


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
    """Computes river water heat potential and temperature profiles for district heating regions."""
    input:
        unpack(input_hera_data),
        onshore_regions=resources("onshore_regions.geojson"),
        dh_areas=resources("dh_areas.geojson"),
    output:
        heat_source_power=resources("heat_source_power_river_water.csv"),
        heat_source_temperature=resources("temp_river_water.nc"),
        heat_source_temperature_temporal_aggregate=resources(
            "temp_river_water_temporal_aggregate.nc"
        ),
        heat_source_energy_temporal_aggregate=resources(
            "heat_source_energy_river_water_temporal_aggregate.nc"
        ),
    log:
        logs("build_river_water_heat_potential.log"),
    benchmark:
        benchmarks("build_river_water_heat_potential")
    threads: 1
    resources:
        mem_mb=20000,
    params:
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        snapshots=config_provider("snapshots"),
        dh_area_buffer=config_provider(
            "sector", "district_heating", "dh_areas", "buffer"
        ),
        enable_heat_source_maps=config_provider("plotting", "enable_heat_source_maps"),
    script:
        scripts(
            "build_surface_water_heat_potentials/build_river_water_heat_potential.py"
        )


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
            + ("_{horizon}" if heat_source_name == "ptes" else "")
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
    """Computes sea water temperature profiles as a district heating heat source per region."""
    input:
        # seawater_temperature=lambda w: input_seawater_temperature(w),
        unpack(input_seawater_temperature),
        onshore_regions=resources("onshore_regions.geojson"),
        dh_areas=resources("dh_areas.geojson"),
    output:
        heat_source_temperature=resources("temp_sea_water.nc"),
        heat_source_temperature_temporal_aggregate=resources(
            "temp_sea_water_temporal_aggregate.nc"
        ),
    log:
        logs("build_sea_water_heat_potential.log"),
    benchmark:
        benchmarks("build_sea_water_heat_potential")
    threads: config["atlite"].get("nprocesses", 4)
    resources:
        mem_mb=10000,
    params:
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        snapshots=config_provider("snapshots"),
        dh_area_buffer=config_provider(
            "sector", "district_heating", "dh_areas", "buffer"
        ),
    script:
        scripts("build_surface_water_heat_potentials/build_sea_water_heat_potential.py")


rule build_cop_profiles:
    """Approximates heat pump coefficient-of-performance profiles for all heat sources and systems."""
    input:
        unpack(input_heat_source_temperature),
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_{horizon}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_{horizon}.nc"
        ),
        temp_soil_total=resources("temp_soil_total.nc"),
        temp_air_total=resources("temp_air_total.nc"),
        temp_ptes_total=resources("ptes_top_temperature_profiles_{horizon}.nc"),
        onshore_regions=resources("onshore_regions.geojson"),
    output:
        cop_profiles=resources("cop_profiles_{horizon}.nc"),
    log:
        logs("build_cop_profiles_{horizon}.log"),
    benchmark:
        benchmarks("build_cop_profiles_{horizon}")
    resources:
        mem_mb=20000,
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
    script:
        scripts("build_cop_profiles/run.py")


rule build_ptes_operations:
    """Builds pit thermal storage top temperature, direct-use and capacity profiles per region."""
    input:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_{horizon}.nc"
        ),
        central_heating_return_temperature_profiles=resources(
            "central_heating_return_temperature_profiles_{horizon}.nc"
        ),
        onshore_regions=resources("onshore_regions.geojson"),
    output:
        ptes_direct_utilisation_profiles=resources(
            "ptes_direct_utilisation_profiles_{horizon}.nc"
        ),
        ptes_top_temperature_profiles=resources(
            "ptes_top_temperature_profiles_{horizon}.nc"
        ),
        ptes_e_max_pu_profiles=resources("ptes_e_max_pu_profiles_{horizon}.nc"),
    log:
        logs("build_ptes_operations_{horizon}.log"),
    benchmark:
        benchmarks("build_ptes_operations_{horizon}")
    resources:
        mem_mb=2000,
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
    script:
        scripts("build_ptes_operations/run.py")


rule build_direct_heat_source_utilisation_profiles:
    """Builds availability profiles for direct heat source use from forward temperature profiles."""
    input:
        central_heating_forward_temperature_profiles=resources(
            "central_heating_forward_temperature_profiles_{horizon}.nc"
        ),
    output:
        direct_heat_source_utilisation_profiles=resources(
            "direct_heat_source_utilisation_profiles_{horizon}.nc"
        ),
    log:
        logs("build_direct_heat_source_utilisation_profiles_{horizon}.log"),
    benchmark:
        benchmarks("build_direct_heat_source_utilisation_profiles_{horizon}")
    resources:
        mem_mb=20000,
    params:
        direct_utilisation_heat_sources=config_provider(
            "sector", "district_heating", "direct_utilisation_heat_sources"
        ),
        limited_heat_sources=config_provider(
            "sector", "district_heating", "limited_heat_sources"
        ),
        snapshots=config_provider("snapshots"),
    script:
        scripts("build_direct_heat_source_utilisation_profiles.py")


rule build_solar_thermal_profiles:
    """Builds solar thermal collector heat generation time series per clustered model region."""
    input:
        pop_layout=resources("pop_layout_total.nc"),
        onshore_regions=resources("onshore_regions.geojson"),
        cutout=lambda w: input_cutout(w, config_provider("solar_thermal", "cutout")(w)),
    output:
        solar_thermal=resources("solar_thermal_total.nc"),
    log:
        logs("build_solar_thermal_profiles_total.log"),
    benchmark:
        benchmarks("build_solar_thermal_profiles/total")
    threads: 16
    resources:
        mem_mb=20000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        solar_thermal=config_provider("solar_thermal"),
    script:
        scripts("build_solar_thermal_profiles.py")


rule build_eurostat_balances:
    """Preprocesses Eurostat energy balances into a tidy table per country, year and carrier."""
    input:
        tsv_gz=rules.retrieve_eurostat_balances.output["tsv_gz"],
    output:
        csv=resources("eurostat_energy_balances.csv"),
    log:
        logs("build_eurostat_balances.log"),
    benchmark:
        benchmarks("build_eurostat_balances")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("build_eurostat_balances.py")


rule build_swiss_energy_balances:
    """Extracts historic Swiss energy balances in TWh per year from the federal spreadsheet."""
    input:
        xlsx=rules.retrieve_swiss_energy_balances.output["xlsx"],
    output:
        csv=resources("switzerland_energy_balances.csv"),
    log:
        logs("build_swiss_energy_balances.log"),
    benchmark:
        benchmarks("build_swiss_energy_balances")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("build_swiss_energy_balances.py")


rule build_co2_totals:
    """Computes historical CO2 emissions per country and sector from EEA and Eurostat data."""
    input:
        co2=rules.retrieve_ghg_emissions.output["csv"],
        eurostat=resources("eurostat_energy_balances.csv"),
    output:
        co2_totals=resources("co2_totals.csv"),
    log:
        logs("build_co2_totals.log"),
    benchmark:
        benchmarks("build_co2_totals")
    threads: 1
    resources:
        mem_mb=1000,
    params:
        countries=config_provider("countries"),
        energy=config_provider("energy"),
        emissions_scope=config_provider("co2_budget", "emissions_scope"),
    script:
        scripts("build_co2_totals.py")


rule build_transformation_output_coke:
    """Extracts coke oven transformation output per country from Eurostat energy balances."""
    input:
        eurostat=resources("eurostat_energy_balances.csv"),
    output:
        transformation_output_coke=resources("transformation_output_coke.csv"),
    log:
        logs("build_transformation_output_coke.log"),
    benchmark:
        benchmarks("build_transformation_output_coke")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("build_transformation_output_coke.py")


rule build_energy_totals:
    """Builds annual energy demand totals per country and sector from JRC IDEES and Eurostat data."""
    input:
        nuts3_shapes=resources("nuts3_shapes.geojson"),
        swiss=resources("switzerland_energy_balances.csv"),
        swiss_transport=lambda w: (
            f"{BFS_ROAD_VEHICLE_STOCK_DATASET['folder']}/vehicle_stock.csv"
            if "CH" in config_provider("countries")(w)
            else []
        ),
        idees=rules.retrieve_jrc_idees.output["directory"],
        district_heat_share="data/district_heat_share.csv",
        eurostat=resources("eurostat_energy_balances.csv"),
        eurostat_households=rules.retrieve_eurostat_household_balances.output["csv"],
    output:
        energy_name=resources("energy_totals.csv"),
        transport_name=resources("transport_data_raw.csv"),
        district_heat_share=resources("district_heat_share.csv"),
        heating_efficiencies=resources("heating_efficiencies.csv"),
    log:
        logs("build_energy_totals.log"),
    benchmark:
        benchmarks("build_energy_totals")
    threads: 16
    resources:
        mem_mb=10000,
    params:
        countries=config_provider("countries"),
        energy=config_provider("energy"),
    script:
        scripts("build_energy_totals.py")


if (COUNTRY_HDD_DATASET := dataset_version("country_hdd"))["source"] in ["build"]:

    # This rule uses one or multiple cutouts.
    # To update the output files to include a new year, e.g. 2025 using an existing cutout,
    # either create a new cutout covering the whole timespan or add another cutout that covers the additional year(s).
    # E.g. cutouts=[<cutout for 1940-2024>, <cutout for 2025-2025>]
    rule build_country_hdd:
        """Computes daily heating degree days per country from ERA5 temperatures for all weather years."""
        input:
            cutouts=["cutouts/europe-1940-2024-era5.nc"],
            country_shapes=resources("country_shapes.geojson"),
        output:
            era5_hdd=f"{COUNTRY_HDD_DATASET['folder']}/era5-HDD-per-country.csv",
        log:
            logs("build_country_hdd.log"),
        benchmark:
            benchmarks("build_country_hdd")
        script:
            scripts("build_country_hdd.py")


rule build_heat_totals:
    """Approximates annual heat demand per country for all weather years via heating degree days."""
    input:
        hdd=f"{COUNTRY_HDD_DATASET['folder']}/era5-HDD-per-country.csv",
        energy_totals=resources("energy_totals.csv"),
    output:
        heat_totals=resources("heat_totals.csv"),
    log:
        logs("build_heat_totals.log"),
    benchmark:
        benchmarks("build_heat_totals")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        scripts("build_heat_totals.py")


rule build_biomass_potentials:
    """Computes biogas and solid biomass potentials per clustered region from JRC ENSPRESO data."""
    input:
        enspreso_biomass=rules.retrieve_enspreso_biomass.output["xlsx"],
        eurostat=resources("eurostat_energy_balances.csv"),
        nuts2=rules.retrieve_eu_nuts_2013.output["shapes_level_2"],
        onshore_regions=resources("onshore_regions.geojson"),
        nuts3_population=ancient(rules.retrieve_nuts3_population.output["gz"]),
        swiss_cantons=lambda w: (
            ancient("data/ch_cantons.csv")
            if "CH" in config_provider("countries")(w)
            else []
        ),
        swiss_population=lambda w: (
            rules.retrieve_bfs_gdp_and_population.output["xlsx"]
            if "CH" in config_provider("countries")(w)
            else []
        ),
        country_shapes=resources("country_shapes.geojson"),
    output:
        biomass_potentials_all=resources("biomass_potentials_all_{horizon}.csv"),
        biomass_potentials=resources("biomass_potentials_{horizon}.csv"),
    log:
        logs("build_biomass_potentials_{horizon}.log"),
    benchmark:
        benchmarks("build_biomass_potentials_{horizon}")
    threads: 8
    resources:
        mem_mb=2000,
    params:
        biomass=config_provider("biomass"),
    script:
        scripts("build_biomass_potentials.py")


rule build_biomass_transport_costs:
    """Converts JRC biomass transport costs per country into EUR per km and MWh."""
    input:
        sc1="data/biomass_transport_costs_supplychain1.csv",
        sc2="data/biomass_transport_costs_supplychain2.csv",
    output:
        biomass_transport_costs=resources("biomass_transport_costs.csv"),
    log:
        logs("build_biomass_transport_costs.log"),
    benchmark:
        benchmarks("build_biomass_transport_costs")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("build_biomass_transport_costs.py")


rule build_co2_sequestration_potentials:
    """Builds geological CO2 sequestration potential shapes from the CO2Stop database."""
    input:
        storage_table=rules.retrieve_co2stop.output["storage_table"],
        storage_map=rules.retrieve_co2stop.output["storage_map"],
        traps_table1=rules.retrieve_co2stop.output["traps_table1"],
        traps_table2=rules.retrieve_co2stop.output["traps_table2"],
        traps_table3=rules.retrieve_co2stop.output["traps_table3"],
        traps_map=rules.retrieve_co2stop.output["traps_map"],
    output:
        resources("co2_sequestration_potentials.geojson"),
    log:
        logs("build_co2_sequestration_potentials.log"),
    benchmark:
        benchmarks("build_co2_sequestration_potentials")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        scripts("build_co2_sequestration_potentials.py")


rule build_clustered_co2_sequestration_potentials:
    """Aggregates geological CO2 sequestration potentials to clustered model regions."""
    input:
        sequestration_potential=resources("co2_sequestration_potentials.geojson"),
        onshore_regions=resources("onshore_regions.geojson"),
        offshore_regions=resources("offshore_regions.geojson"),
    output:
        sequestration_potential=resources("co2_sequestration_potential.csv"),
    log:
        logs("build_clustered_co2_sequestration_potentials.log"),
    benchmark:
        benchmarks("build_clustered_co2_sequestration_potentials")
    threads: 1
    resources:
        mem_mb=4000,
    params:
        sequestration_potential=config_provider(
            "sector", "regional_co2_sequestration_potential"
        ),
    script:
        scripts("build_clustered_co2_sequestration_potentials.py")


rule build_salt_cavern_potentials:
    """Builds hydrogen storage potentials in salt caverns per region split by onshore and offshore."""
    input:
        salt_caverns=rules.retrieve_h2_salt_caverns.output["geojson"],
        onshore_regions=resources("onshore_regions.geojson"),
        offshore_regions=resources("offshore_regions.geojson"),
    output:
        h2_cavern_potential=resources("salt_cavern_potentials.csv"),
    log:
        logs("build_salt_cavern_potentials.log"),
    benchmark:
        benchmarks("build_salt_cavern_potentials")
    threads: 1
    resources:
        mem_mb=2000,
    script:
        scripts("build_salt_cavern_potentials.py")


rule build_ammonia_production:
    """Extracts historical annual ammonia production per country from USGS statistics."""
    input:
        usgs=rules.retrieve_nitrogen_statistics.output["xlsx"],
    output:
        ammonia_production=resources("ammonia_production.csv"),
    log:
        logs("build_ammonia_production.log"),
    benchmark:
        benchmarks("build_ammonia_production")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("build_ammonia_production.py")


rule build_industry_sector_ratios:
    """Builds best-case specific energy consumption per carrier and industry from JRC IDEES."""
    input:
        ammonia_production=resources("ammonia_production.csv"),
        idees=rules.retrieve_jrc_idees.output["directory"],
    output:
        industry_sector_ratios=resources("industry_sector_ratios.csv"),
    log:
        logs("build_industry_sector_ratios.log"),
    benchmark:
        benchmarks("build_industry_sector_ratios")
    threads: 1
    resources:
        mem_mb=1000,
    params:
        industry=config_provider("industry"),
        ammonia=config_provider("sector", "ammonia", default=False),
    script:
        scripts("build_industry_sector_ratios.py")


rule build_industry_sector_ratios_intermediate:
    """Interpolates specific industrial energy consumption between today and best-in-class."""
    input:
        industry_sector_ratios=resources("industry_sector_ratios.csv"),
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    output:
        industry_sector_ratios=resources("industry_sector_ratios_{horizon}.csv"),
    log:
        logs("build_industry_sector_ratios_{horizon}.log"),
    benchmark:
        benchmarks("build_industry_sector_ratios_{horizon}")
    threads: 1
    resources:
        mem_mb=1000,
    params:
        industry=config_provider("industry"),
    script:
        scripts("build_industry_sector_ratios_intermediate.py")


rule build_industrial_production_per_country:
    """Builds historical industrial production per country from JRC IDEES and Eurostat data."""
    input:
        ch_industrial_production="data/ch_industrial_production_per_subsector.csv",
        ammonia_production=resources("ammonia_production.csv"),
        eurostat=resources("eurostat_energy_balances.csv"),
        jrc=rules.retrieve_jrc_idees.output["directory"],
    output:
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    log:
        logs("build_industrial_production_per_country.log"),
    benchmark:
        benchmarks("build_industrial_production_per_country")
    threads: 8
    resources:
        mem_mb=2000,
    params:
        industry=config_provider("industry"),
        countries=config_provider("countries"),
    script:
        scripts("build_industrial_production_per_country.py")


rule build_industrial_production_per_country_tomorrow:
    """Projects future industrial production per country from recycling and primary shares."""
    input:
        industrial_production_per_country=resources(
            "industrial_production_per_country.csv"
        ),
    output:
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{horizon}.csv"
        ),
    log:
        logs("build_industrial_production_per_country_tomorrow_{horizon}.log"),
    benchmark:
        (benchmarks("build_industrial_production_per_country_tomorrow_{horizon}"))
    threads: 1
    resources:
        mem_mb=1000,
    params:
        industry=config_provider("industry"),
    script:
        scripts("build_industrial_production_per_country_tomorrow.py")


rule build_industrial_distribution_key:
    """Builds nodal distribution keys per industry sector from Hotmaps industrial sites."""
    input:
        onshore_regions=resources("onshore_regions.geojson"),
        clustered_pop_layout=resources("pop_layout.csv"),
        hotmaps=rules.retrieve_hotmaps_industrial_sites.output["csv"],
        gem_gspt=rules.retrieve_gem_steel_plant_tracker.output["xlsx"],
        gem_gcpt=rules.retrieve_gem_cement_concrete_tracker.output["xlsx"],
        ammonia="data/ammonia_plants.csv",
        refineries_supplement="data/refineries-noneu.csv",
    output:
        industrial_distribution_key=resources("industrial_distribution_key.csv"),
    log:
        logs("build_industrial_distribution_key.log"),
    benchmark:
        benchmarks("build_industrial_distribution_key")
    threads: 1
    resources:
        mem_mb=1000,
    params:
        hotmaps_locate_missing=config_provider(
            "industry", "hotmaps_locate_missing", default=False
        ),
        countries=config_provider("countries"),
    script:
        scripts("build_industrial_distribution_key.py")


rule build_industrial_production_per_node:
    """Distributes industrial production per country to model regions with distribution keys."""
    input:
        industrial_distribution_key=resources("industrial_distribution_key.csv"),
        industrial_production_per_country_tomorrow=resources(
            "industrial_production_per_country_tomorrow_{horizon}.csv"
        ),
    output:
        industrial_production_per_node=resources("industrial_production_{horizon}.csv"),
    log:
        logs("build_industrial_production_per_node_{horizon}.log"),
    benchmark:
        (benchmarks("build_industrial_production_per_node_{horizon}"))
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("build_industrial_production_per_node.py")


rule build_industrial_energy_demand_per_node:
    """Computes industrial energy demand per carrier and model region from production and ratios."""
    input:
        industry_sector_ratios=resources("industry_sector_ratios_{horizon}.csv"),
        industrial_production_per_node=resources("industrial_production_{horizon}.csv"),
        industrial_energy_demand_per_node_today=resources(
            "industrial_energy_demand_today.csv"
        ),
    output:
        industrial_energy_demand_per_node=resources(
            "industrial_energy_demand_{horizon}.csv"
        ),
    log:
        logs("build_industrial_energy_demand_per_node_{horizon}.log"),
    benchmark:
        (benchmarks("build_industrial_energy_demand_per_node_{horizon}"))
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("build_industrial_energy_demand_per_node.py")


rule build_industrial_energy_demand_per_country_today:
    """Computes today's industrial energy demand per country and sector from JRC IDEES."""
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
    log:
        logs("build_industrial_energy_demand_per_country_today.log"),
    benchmark:
        benchmarks("build_industrial_energy_demand_per_country_today")
    threads: 8
    resources:
        mem_mb=2000,
    params:
        countries=config_provider("countries"),
        industry=config_provider("industry"),
        ammonia=config_provider("sector", "ammonia", default=False),
    script:
        scripts("build_industrial_energy_demand_per_country_today.py")


rule build_industrial_energy_demand_per_node_today:
    """Distributes today's industrial energy demand per country to model regions."""
    input:
        industrial_distribution_key=resources("industrial_distribution_key.csv"),
        industrial_energy_demand_per_country_today=resources(
            "industrial_energy_demand_per_country_today.csv"
        ),
    output:
        industrial_energy_demand_per_node_today=resources(
            "industrial_energy_demand_today.csv"
        ),
    log:
        logs("build_industrial_energy_demand_per_node_today.log"),
    benchmark:
        benchmarks("build_industrial_energy_demand_per_node_today")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        scripts("build_industrial_energy_demand_per_node_today.py")


rule build_retro_cost:
    """Computes building retrofit costs and space heating savings per region and building type."""
    input:
        building_stock="data/retro/data_building_stock.csv",
        data_tabula=rules.retrieve_tabula_calculator.output["xlsx"],
        air_temperature=resources("temp_air_total.nc"),
        u_values_PL="data/retro/u_values_poland.csv",
        tax_w="data/retro/electricity_taxes_eu.csv",
        construction_index="data/retro/comparative_level_investment.csv",
        floor_area_missing="data/retro/floor_area_missing.csv",
        clustered_pop_layout=resources("pop_layout.csv"),
        cost_germany="data/retro/retro_cost_germany.csv",
        window_assumptions="data/retro/window_assumptions.csv",
    output:
        retro_cost=resources("retro_cost.csv"),
        floor_area=resources("floor_area.csv"),
    log:
        logs("build_retro_cost.log"),
    benchmark:
        benchmarks("build_retro_cost")
    resources:
        mem_mb=1000,
    params:
        retrofitting=config_provider("sector", "retrofitting"),
        countries=config_provider("countries"),
    script:
        scripts("build_retro_cost.py")


rule build_population_weighted_energy_totals:
    """Distributes country-level energy demand totals to model regions by population."""
    input:
        energy_totals=resources("{kind}_totals.csv"),
        clustered_pop_layout=resources("pop_layout.csv"),
    output:
        resources("pop_weighted_{kind}_totals.csv"),
    log:
        logs("build_population_weighted_{kind}_totals.log"),
    benchmark:
        benchmarks("build_population_weighted_{kind}_totals")
    threads: 1
    resources:
        mem_mb=2000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    script:
        scripts("build_population_weighted_energy_totals.py")


rule build_shipping_demand:
    """Builds regional international shipping energy demand from port outflow volumes."""
    input:
        ports=rules.retrieve_attributed_ports.output["json"],
        scope=resources("europe_shape.geojson"),
        regions=resources("onshore_regions.geojson"),
        demand=resources("energy_totals.csv"),
    output:
        resources("shipping_demand.csv"),
    log:
        logs("build_shipping_demand.log"),
    benchmark:
        benchmarks("build_shipping_demand")
    threads: 1
    resources:
        mem_mb=2000,
    params:
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    script:
        scripts("build_shipping_demand.py")


if MOBILITY_PROFILES_DATASET["source"] in ["build"]:

    rule build_mobility_profiles:
        """Builds weekly road transport profiles from German BASt vehicle count data."""
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
        log:
            logs("build_mobility_profiles.log"),
        benchmark:
            benchmarks("build_mobility_profiles")
        threads: 1
        resources:
            mem_mb=5000,
        script:
            scripts("build_mobility_profiles.py")


rule build_transport_demand:
    """Builds land transport demand and electric vehicle availability profiles per region."""
    input:
        network=resources("networks/clustered.nc"),
        clustered_pop_layout=resources("pop_layout.csv"),
        pop_weighted_energy_totals=resources("pop_weighted_energy_totals.csv"),
        transport_data_raw=resources("transport_data_raw.csv"),
        traffic_data_KFZ=f"{MOBILITY_PROFILES_DATASET['folder']}/kfz.csv",
        traffic_data_Pkw=f"{MOBILITY_PROFILES_DATASET['folder']}/pkw.csv",
        temp_air_total=resources("temp_air_total.nc"),
    output:
        transport_demand=resources("transport_demand.csv"),
        transport_data=resources("transport_data.csv"),
        avail_profile=resources("avail_profile.csv"),
        dsm_profile=resources("dsm_profile.csv"),
    log:
        logs("build_transport_demand.log"),
    benchmark:
        benchmarks("build_transport_demand")
    threads: 1
    resources:
        mem_mb=2000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        sector=config_provider("sector"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    script:
        scripts("build_transport_demand.py")


rule build_district_heat_share:
    """Builds district heating shares per region and investment year from urban population."""
    input:
        district_heat_share=resources("district_heat_share.csv"),
        clustered_pop_layout=resources("pop_layout.csv"),
    output:
        district_heat_share=resources("district_heat_share_{horizon}.csv"),
    log:
        logs("build_district_heat_share_{horizon}.log"),
    benchmark:
        benchmarks("build_district_heat_share_{horizon}")
    threads: 1
    resources:
        mem_mb=1000,
    params:
        sector=config_provider("sector"),
        energy_totals_year=config_provider("energy", "energy_totals_year"),
    script:
        scripts("build_district_heat_share.py")


rule build_existing_heating_distribution:
    """Distributes existing heating capacities to regions and sectors by population."""
    input:
        existing_heating="data/existing_infrastructure/existing_heating_raw.csv",
        clustered_pop_layout=resources("pop_layout.csv"),
        clustered_pop_energy_layout=resources("pop_weighted_energy_totals.csv"),
        district_heat_share=resources("district_heat_share_{horizon}.csv"),
    output:
        existing_heating_distribution=resources(
            "existing_heating_distribution_{horizon}.csv"
        ),
    log:
        logs("build_existing_heating_distribution_{horizon}.log"),
    benchmark:
        benchmarks("build_existing_heating_distribution_{horizon}")
    threads: 1
    resources:
        mem_mb=2000,
    params:
        baseyear=config_provider("planning_horizons", default=0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
    script:
        scripts("build_existing_heating_distribution.py")


rule time_aggregation:
    """Computes snapshot weightings for the time aggregation of the sector-coupled network."""
    input:
        network=resources("networks/clustered.nc"),
        electricity_demand=resources("electricity_demand.nc"),
        profiles=lambda w: [
            resources(f"profile_{tech}.nc")
            for tech in config_provider("electricity", "renewable_carriers")(w)
            if tech != "hydro"
        ],
        hydro_profile=lambda w: (
            resources("profile_hydro.nc")
            if "hydro" in config_provider("electricity", "renewable_carriers")(w)
            else []
        ),
        hourly_heat_demand_total=lambda w: (
            resources("hourly_heat_demand_total.nc")
            if config_provider("sector", "enabled")(w)
            and config_provider("sector", "heating")(w)
            else []
        ),
        solar_thermal_total=lambda w: (
            resources("solar_thermal_total.nc")
            if config_provider("sector", "enabled")(w)
            and config_provider("sector", "solar_thermal")(w)
            else []
        ),
        # TODO: add cop and transport profiles, search for others in prepare_sectcor (like master considers them) 
    output:
        snapshot_weightings=resources("snapshot_weightings.csv"),
    log:
        logs("time_aggregation_elec.log"),
    benchmark:
        benchmarks("time_aggregation")
    threads: 1
    resources:
        mem_mb=5000,
    params:
        time_resolution=config_provider("clustering", "temporal"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        solver_name=config_provider("solving", "solver", "name"),
    script:
        scripts("time_aggregation.py")


def input_profile_offwind(w):
    return {
        f"profile_{tech}": resources("profile_" + tech + ".nc")
        for tech in ["offwind-ac", "offwind-dc", "offwind-float"]
        if (tech in config_provider("electricity", "renewable_carriers")(w))
    }


rule build_egs_potentials:
    """Builds enhanced geothermal capacity potentials and costs per region from gridded data."""
    input:
        egs_cost="data/egs_costs.json",
        regions=resources("onshore_regions.geojson"),
        air_temperature=(
            resources("temp_air_total.nc")
            if config_provider("sector", "enhanced_geothermal", "var_cf")
            else []
        ),
    output:
        egs_potentials=resources("egs_potentials.csv"),
        egs_overlap=resources("egs_overlap.csv"),
        egs_capacity_factors=resources("egs_capacity_factors.csv"),
    log:
        logs("build_egs_potentials.log"),
    benchmark:
        benchmarks("build_egs_potentials")
    threads: 1
    resources:
        mem_mb=2000,
    params:
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        sector=config_provider("sector"),
        costs=config_provider("costs"),
    script:
        scripts("build_egs_potentials.py")


def input_heat_source_power(w):

    return {
        heat_source_name: resources("heat_source_power_" + heat_source_name + ".csv")
        for heat_source_name in config_provider(
            "sector", "heat_pump_sources", "urban central"
        )(w)
        if heat_source_name
        in config_provider("sector", "district_heating", "limited_heat_sources")(
            w
        ).keys()
    }
