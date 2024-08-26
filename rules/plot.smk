# SPDX-FileCopyrightText: 2023- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule plot_power_network_unclustered:
    input:
        network=resources("networks/elec.nc"),
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/power-network-unclustered"), ".png", ".pdf"),
    script:
        "../scripts/plot_power_network_unclustered.py"


rule plot_gas_network_unclustered:
    input:
        gas_network=resources("gas_network.csv"),
        gem="data/gem/Europe-Gas-Tracker-2024-05.xlsx",
        entry="data/gas_network/scigrid-gas/data/IGGIELGN_BorderPoints.geojson",
        storage="data/gas_network/scigrid-gas/data/IGGIELGN_Storages.geojson",
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/gas-network-unclustered"), ".png", ".pdf"),
    script:
        "../scripts/plot_gas_network_unclustered.py"


rule plot_renewable_potential_unclustered:
    input:
        **{
            f"profile_{tech}": resources(f"profile_{tech}.nc")
            for tech in config["electricity"]["renewable_carriers"]
        },
        regions_onshore=resources("regions_onshore.geojson"),
        regions_offshore=resources("regions_offshore.geojson"),
        rc="matplotlibrc",
    output:
        wind=multiext(resources("graphics/wind-energy-density"), ".png", ".pdf"),
        solar=multiext(resources("graphics/solar-energy-density"), ".png", ".pdf"),
        wind_cf=multiext(resources("graphics/wind-capacity-factor"), ".png", ".pdf"),
        solar_cf=multiext(resources("graphics/solar-capacity-factor"), ".png", ".pdf"),
    script:
        "../scripts/plot_renewable_potential_unclustered.py"


rule plot_weather_data_map:
    input:
        cutout=f"cutouts/" + CDIR + config["atlite"]["default_cutout"] + ".nc",
        rc="matplotlibrc",
    output:
        irradiation=multiext(
            resources("graphics/weather-map-irradiation"), ".png", ".pdf"
        ),
        runoff=multiext(resources("graphics/weather-map-runoff"), ".png", ".pdf"),
        temperature=multiext(
            resources("graphics/weather-map-temperature"), ".png", ".pdf"
        ),
        wind=multiext(resources("graphics/weather-map-wind"), ".png", ".pdf"),
    script:
        "../scripts/plot_weather_data_map.py"


rule plot_industrial_sites:
    input:
        hotmaps="data/Industrial_Database.csv",
        countries=resources("country_shapes.geojson"),
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/industrial-sites"), ".png", ".pdf"),
    script:
        "../scripts/plot_industrial_sites.py"


rule plot_powerplants:
    input:
        powerplants=resources("powerplants.csv"),
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/powerplants"), ".png", ".pdf"),
    script:
        "../scripts/plot_powerplants.py"


rule plot_salt_caverns_unclustered:
    input:
        caverns="data/h2_salt_caverns_GWh_per_sqkm.geojson",
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/salt-caverns"), ".png", ".pdf"),
    script:
        "../scripts/plot_salt_caverns_unclustered.py"


rule plot_salt_caverns_clustered:
    input:
        caverns=resources("salt_cavern_potentials_s{simpl}_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        onshore=multiext(
            resources("graphics/salt-caverns-s{simpl}-{clusters}-onshore"),
            ".png",
            ".pdf",
        ),
        nearshore=multiext(
            resources("graphics/salt-caverns-s{simpl}-{clusters}-nearshore"),
            ".png",
            ".pdf",
        ),
        offshore=multiext(
            resources("graphics/salt-caverns-s{simpl}-{clusters}-offshore"),
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_salt_caverns_clustered.py"


rule plot_biomass_potentials:
    input:
        biomass=resources("biomass_potentials_s{simpl}_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        solid_biomass=multiext(
            resources("graphics/biomass-potentials-s{simpl}-{clusters}-solid_biomass"),
            ".png",
            ".pdf",
        ),
        not_included=multiext(
            resources("graphics/biomass-potentials-s{simpl}-{clusters}-not_included"),
            ".png",
            ".pdf",
        ),
        biogas=multiext(
            resources("graphics/biomass-potentials-s{simpl}-{clusters}-biogas"),
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_biomass_potentials.py"


rule plot_choropleth_capacity_factors:
    input:
        network=resources("networks/elec_s{simpl}_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(resources("graphics/capacity-factor/s{simpl}-{clusters}")),
    script:
        "../scripts/plot_choropleth_capacity_factors.py"


rule plot_choropleth_capacity_factors_sector:
    input:
        cop_soil=resources("cop_soil_total_elec_s{simpl}_{clusters}.nc"),
        cop_air=resources("cop_air_total_elec_s{simpl}_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(resources("graphics/capacity-factor-sector/s{simpl}-{clusters}")),
    script:
        "../scripts/plot_choropleth_capacity_factors_sector.py"


rule plot_choropleth_capacities:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/p_nom_opt/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_capacities.py"


rule plot_choropleth_prices:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/prices/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_prices.py"


rule plot_choropleth_potential_used:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/potential_used/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_potential_used.py"


rule plot_choropleth_demand:
    input:
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        industrial_demand=resources("industrial_energy_demand_elec_s{simpl}_{clusters}_{planning_horizons}.csv"),
        shipping_demand=resources("shipping_demand_s{simpl}_{clusters}.csv"),
        nodal_energy_totals=resources("pop_weighted_energy_totals_s{simpl}_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/regional_demand/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_demand.py"


rule plot_balance_timeseries:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    threads: 12
    output:
        directory(
            RESULTS
            + "graphics/balance_timeseries/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_balance_timeseries.py"


rule plot_heatmap_timeseries:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    threads: 12
    output:
        directory(
            RESULTS
            + "graphics/heatmap_timeseries/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_heatmap_timeseries.py"


rule plot_heatmap_timeseries_resources:
    input:
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    threads: 12
    output:
        directory(
            RESULTS
            + "graphics/heatmap_timeseries_resources/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_heatmap_timeseries_resources.py"


rule plot_import_options:
    input:
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        entrypoints=resources("gas_input_locations_s{simpl}_{clusters}_simplified.csv"),
        imports="data/imports/results.csv",
        rc="matplotlibrc",
    output:
        map=multiext(
            RESULTS
            + "graphics/import_options_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
            ".png",
            ".pdf",
        ),
        distribution=multiext(
            RESULTS
            + "graphics/import_options_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}-distribution",
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_import_options.py"


rule plot_import_world_map:
    input:
        imports="data/imports/results.csv",
        profiles="data/imports/combined_weighted_generator_timeseries.nc",
        gadm_arg=storage(
            "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_ARG.gpkg",
            keep_local=True,
        ),
        copernicus_glc=storage(
            "https://zenodo.org/records/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
            keep_local=True,
        ),
        wdpa="data/WDPA.gpkg",
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/import_world_map"), ".png", ".pdf"),
    script:
        "../scripts/plot_import_world_map.py"


rule plot_import_networks:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_elec_s{simpl}_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        multiext(
            RESULTS
            + "graphics/import_networks/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_import_network.py"


rule plot_import_shares:
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    output:
        multiext(
            RESULTS
            + "graphics/import_shares/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_import_shares.py"


rule plot_all_resources:
    input:
        rules.plot_power_network_unclustered.output,
        rules.plot_gas_network_unclustered.output,
        rules.plot_renewable_potential_unclustered.output,
        rules.plot_weather_data_map.output,
        rules.plot_industrial_sites.output,
        rules.plot_powerplants.output,
        rules.plot_salt_caverns_unclustered.output,
        rules.plot_import_world_map.output,
        expand(
            resources("graphics/salt-caverns-s{simpl}-{clusters}-onshore.pdf"),
            **config["scenario"],
        ),
        expand(
            resources("graphics/power-network-{clusters}.pdf"), **config["scenario"]
        ),
        expand(
            resources("graphics/biomass-potentials-s{simpl}-{clusters}-biogas.pdf"),
            **config["scenario"],
        ),
        expand(
            resources("graphics/capacity-factor/s{simpl}-{clusters}"),
            **config["scenario"],
        ),
        expand(
            resources("graphics/capacity-factor-sector/s{simpl}-{clusters}"),
            **config["scenario"],
        ),


rule plot_all_results_single:
    input:
        RESULTS
        + "graphics/p_nom_opt/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/prices/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/potential_used/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/balance_timeseries/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/heatmap_timeseries/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/heatmap_timeseries_resources/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/regional_demand/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
        RESULTS
        + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf",
        RESULTS
        + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf",
        RESULTS
        + "graphics/import_options_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
        RESULTS
        + "graphics/import_networks/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
        RESULTS
        + "graphics/import_shares/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
    output:
        RESULTS
        + "graphics/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.touch",
    shell:
        "touch {output}"


rule plot_all_results:
    input:
        expand(
            RESULTS
            + "graphics/s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.touch",
            **config["scenario"],
        ),
