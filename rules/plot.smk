# SPDX-FileCopyrightText: 2023- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule plot_power_network_unclustered:
    input:
        network=resources("networks/base.nc"),
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
        caverns=resources("salt_cavern_potentials_s_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        onshore=multiext(
            resources("graphics/salt-caverns-s-{clusters}-onshore"),
            ".png",
            ".pdf",
        ),
        nearshore=multiext(
            resources("graphics/salt-caverns-s-{clusters}-nearshore"),
            ".png",
            ".pdf",
        ),
        offshore=multiext(
            resources("graphics/salt-caverns-s-{clusters}-offshore"),
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_salt_caverns_clustered.py"


rule plot_biomass_potentials:
    input:
        biomass=resources("biomass_potentials_s_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        solid_biomass=multiext(
            resources("graphics/biomass-potentials-s-{clusters}-solid_biomass"),
            ".png",
            ".pdf",
        ),
        not_included=multiext(
            resources("graphics/biomass-potentials-s-{clusters}-not_included"),
            ".png",
            ".pdf",
        ),
        biogas=multiext(
            resources("graphics/biomass-potentials-s-{clusters}-biogas"),
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_biomass_potentials.py"


rule plot_choropleth_capacity_factors:
    input:
        network=resources("networks/base_s_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(resources("graphics/capacity-factor/s-{clusters}")),
    script:
        "../scripts/plot_choropleth_capacity_factors.py"


rule plot_choropleth_capacity_factors_sector:
    input:
        cop_soil=resources("cop_soil_total_base_s_{clusters}.nc"),
        cop_air=resources("cop_air_total_base_s_{clusters}.nc"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(resources("graphics/capacity-factor-sector/s-{clusters}")),
    script:
        "../scripts/plot_choropleth_capacity_factors_sector.py"


rule plot_choropleth_capacities:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/p_nom_opt/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_capacities.py"

rule plot_dh_share:
    input:
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        dh_share=resources("district_heat_share_base_s_{clusters}_{planning_horizons}.csv"),
        rc="matplotlibrc",
    output:
        multiext(
            resources("graphics/dh-share-s-{clusters}-{planning_horizons}"),
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_dh_share.py"

rule plot_choropleth_prices:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/prices/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_prices.py"


rule plot_choropleth_potential_used:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        regions_offshore=resources("regions_offshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/potential_used/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_potential_used.py"


rule plot_choropleth_demand:
    input:
        network=RESULTS
        + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        industrial_demand=resources(
            "industrial_energy_demand_base_s_{clusters}_{planning_horizons}.csv"
        ),
        shipping_demand=resources("shipping_demand_s_{clusters}.csv"),
        nodal_energy_totals=resources("pop_weighted_energy_totals_s_{clusters}.csv"),
        regions_onshore=resources("regions_onshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        directory(
            RESULTS
            + "graphics/regional_demand/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_choropleth_demand.py"


rule plot_backup_power_map:
    params:
        plotting=config_provider("plotting"),
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        RESULTS + 
        "graphics/backup_map/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
    script:
        "../scripts/plot_backup_power_map.py"


rule plot_balance_timeseries:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    threads: 12
    output:
        directory(
            RESULTS
            + "graphics/balance_timeseries/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_balance_timeseries.py"


rule plot_heatmap_timeseries:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    threads: 12
    output:
        directory(
            RESULTS
            + "graphics/heatmap_timeseries/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_heatmap_timeseries.py"


rule plot_heatmap_timeseries_resources:
    input:
        network=RESULTS
        + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    threads: 12
    output:
        directory(
            RESULTS
            + "graphics/heatmap_timeseries_resources/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        ),
    script:
        "../scripts/plot_heatmap_timeseries_resources.py"


rule plot_import_options:
    input:
        network=RESULTS
        + "prenetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        entrypoints=resources("gas_input_locations_s_{clusters}_simplified.csv"),
        imports="data/imports/results.parquet",
        rc="matplotlibrc",
    output:
        map=multiext(
            RESULTS
            + "graphics/import_options_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
            ".png",
            ".pdf",
        ),
        distribution=multiext(
            RESULTS
            + "graphics/import_options_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}-distribution",
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_import_options.py"


rule retrieve_gadm_argentina:
    input:
        storage("https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_ARG.gpkg"),
    output:
        "data/gadm/gadm41_ARG.gpkg",
    run:
        move(input[0], output[0])


rule retrieve_naturalearth_lowres_countries:
    input:
        storage(
            "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip"
        ),
    params:
        zip="data/naturalearth/ne_110m_admin_0_countries.zip",
    output:
        countries="data/naturalearth/ne_110m_admin_0_countries.shp",
    run:
        move(input[0], params["zip"])
        output_folder = Path(output["countries"]).parent
        unpack_archive(params["zip"], output_folder)
        os.remove(params["zip"])


rule plot_import_world_map:
    input:
        imports="data/imports/results.parquet",
        profiles="data/imports/combined_weighted_generator_timeseries.nc",
        gadm_arg="data/gadm/gadm41_ARG.gpkg",
        countries="data/naturalearth/ne_110m_admin_0_countries.shp",
        copernicus_glc="data/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        wdpa="data/WDPA.gpkg",
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/import_world_map"), ".png", ".pdf"),
    script:
        "../scripts/plot_import_world_map.py"


rule plot_import_world_map_hydrogen:
    input:
        imports="data/imports/results.parquet",
        countries="data/naturalearth/ne_110m_admin_0_countries.shp",
        rc="matplotlibrc",
    output:
        multiext(resources("graphics/import_world_map_hydrogen"), ".png", ".pdf"),
    script:
        "../scripts/plot_import_world_map_hydrogen.py"


rule plot_import_networks:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=resources("regions_onshore_base_s_{clusters}.geojson"),
        rc="matplotlibrc",
    output:
        multiext(
            RESULTS
            + "graphics/import_networks/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_import_network.py"


rule plot_import_shares:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        rc="matplotlibrc",
    output:
        multiext(
            RESULTS
            + "graphics/import_shares/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
            ".png",
            ".pdf",
        ),
    script:
        "../scripts/plot_import_shares.py"


rule plot_import_sankey:
    input:
        network=RESULTS
        + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    output:
        RESULTS
        + "graphics/import_sankey/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
    script:
        "../scripts/plot_import_sankey.py"


rule plot_import_supply_curve:
    input:
        imports="data/imports/results.parquet",
        rc="matplotlibrc",
    output:
        hvdc_to_elec=multiext(resources("graphics/import_supply_curve_hvdc-to-elec"), ".png", ".pdf"),
        pipeline_h2=multiext(resources("graphics/import_supply_curve_pipeline-h2"), ".png", ".pdf"),
        shipping_lh2=multiext(resources("graphics/import_supply_curve_shipping-lh2"), ".png", ".pdf"),
        shipping_lnh3=multiext(resources("graphics/import_supply_curve_shipping-lnh3"), ".png", ".pdf"),
        shipping_lch4=multiext(resources("graphics/import_supply_curve_shipping-lch4"), ".png", ".pdf"),
        shipping_meoh=multiext(resources("graphics/import_supply_curve_shipping-meoh"), ".png", ".pdf"),
        shipping_ftfuel=multiext(resources("graphics/import_supply_curve_shipping-ftfuel"), ".png", ".pdf"),
        shipping_hbi=multiext(resources("graphics/import_supply_curve_shipping-hbi"), ".png", ".pdf"),
        shipping_steel=multiext(resources("graphics/import_supply_curve_shipping-steel"), ".png", ".pdf"),
    script:
        "../scripts/plot_import_supply_curve.py"

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
            resources("graphics/salt-caverns-s-{clusters}-onshore.pdf"),
            **config["scenario"],
        ),
        expand(resources("graphics/power-network-{clusters}.pdf"), **config["scenario"]),
        expand(
            resources("graphics/biomass-potentials-s-{clusters}-biogas.pdf"),
            **config["scenario"],
        ),
        expand(
            resources("graphics/capacity-factor/s-{clusters}"),
            **config["scenario"],
        ),
        expand(
            resources("graphics/capacity-factor-sector/s-{clusters}"),
            **config["scenario"],
        ),


rule plot_all_results_single:
    input:
        RESULTS
        + "graphics/p_nom_opt/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/prices/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/potential_used/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/balance_timeseries/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/heatmap_timeseries/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/heatmap_timeseries_resources/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "graphics/regional_demand/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}",
        RESULTS
        + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
        RESULTS
        + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf",
        RESULTS
        + "maps/base_s_{clusters}_l{ll}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf",
        RESULTS
        + "graphics/import_options_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
        RESULTS
        + "graphics/import_networks/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
        RESULTS
        + "graphics/import_shares/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.pdf",
    output:
        RESULTS
        + "graphics/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.touch",
    shell:
        "touch {output}"


rule plot_all_results:
    input:
        expand(
            RESULTS
            + "graphics/s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.touch",
            **config["scenario"],
        ),
