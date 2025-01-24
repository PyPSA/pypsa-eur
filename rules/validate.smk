# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

PRODUCTION_PLOTS = [
    "production_bar",
    "production_deviation_bar",
    "seasonal_operation_area",
]
CROSS_BORDER_PLOTS = ["trade_time_series", "cross_border_bar"]
PRICES_PLOTS = ["price_bar", "price_line"]


rule build_electricity_production:
    """
    This rule builds the electricity production for each country and technology from ENTSO-E data.
    The data is used for validation of the optimization results.
    """
    params:
        snapshots=config_provider("snapshots"),
        countries=config_provider("countries"),
    output:
        resources("historical_electricity_production.csv"),
    log:
        logs("build_electricity_production.log"),
    resources:
        mem_mb=5000,
    script:
        "../scripts/build_electricity_production.py"


rule build_cross_border_flows:
    """
    This rule builds the cross-border flows from ENTSO-E data.
    The data is used for validation of the optimization results.
    """
    params:
        snapshots=config_provider("snapshots"),
        countries=config_provider("countries"),
    input:
        network=resources("networks/base.nc"),
    output:
        resources("historical_cross_border_flows.csv"),
    log:
        logs("build_cross_border_flows.log"),
    resources:
        mem_mb=5000,
    script:
        "../scripts/build_cross_border_flows.py"


rule build_electricity_prices:
    """
    This rule builds the electricity prices from ENTSO-E data.
    The data is used for validation of the optimization results.
    """
    params:
        snapshots=config_provider("snapshots"),
        countries=config_provider("countries"),
    output:
        resources("historical_electricity_prices.csv"),
    log:
        logs("build_electricity_prices.log"),
    resources:
        mem_mb=5000,
    script:
        "../scripts/build_electricity_prices.py"


rule plot_validation_electricity_production:
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
        electricity_production=resources("historical_electricity_production.csv"),
    output:
        **{
            plot: RESULTS
            + f"figures/validation_{plot}_base_s_{{clusters}}_elec_{{opts}}.pdf"
            for plot in PRODUCTION_PLOTS
        },
        plots_touch=RESULTS
        + "figures/.validation_production_plots_base_s_{clusters}_elec_{opts}",
    script:
        "../scripts/plot_validation_electricity_production.py"


rule plot_validation_cross_border_flows:
    params:
        countries=config_provider("countries"),
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
        cross_border_flows=resources("historical_cross_border_flows.csv"),
    output:
        **{
            plot: RESULTS
            + f"figures/validation_{plot}_base_s_{{clusters}}_elec_{{opts}}.pdf"
            for plot in CROSS_BORDER_PLOTS
        },
        plots_touch=RESULTS
        + "figures/.validation_cross_border_plots_base_s_{clusters}_elec_{opts}",
    script:
        "../scripts/plot_validation_cross_border_flows.py"


rule plot_validation_electricity_prices:
    input:
        network=RESULTS + "networks/base_s_{clusters}_elec_{opts}.nc",
        electricity_prices=resources("historical_electricity_prices.csv"),
    output:
        **{
            plot: RESULTS
            + f"figures/validation_{plot}_base_s_{{clusters}}_elec_{{opts}}.pdf"
            for plot in PRICES_PLOTS
        },
        plots_touch=RESULTS
        + "figures/.validation_prices_plots_base_s_{clusters}_elec_{opts}",
    script:
        "../scripts/plot_validation_electricity_prices.py"
