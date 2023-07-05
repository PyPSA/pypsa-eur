# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


rule build_electricity_production:
    """
    This rule builds the electricity production for each country and technology from ENTSO-E data.
    The data is used for validation of the optimization results.
    """
    params:
        snapshots=config["snapshots"],
        countries=config["countries"],
    output:
        RESOURCES + "historical_electricity_production.csv",
    log:
        LOGS + "build_electricity_production.log",
    resources:
        mem_mb=5000,
    script:
        "../scripts/build_electricity_production.py"


PLOTS = ["production_bar", "production_deviation_bar", "seasonal_operation_area"]


rule plot_electricity_production:
    input:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        electricity_production=RESOURCES + "historical_electricity_production.csv",
    output:
        **{
            plot: RESULTS
            + f"figures/validation_{plot}_elec_s{{simpl}}_{{clusters}}_ec_l{{ll}}_{{opts}}.pdf"
            for plot in PLOTS
        },
        plots_touch=RESULTS
        + "figures/.validation_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
    script:
        "../scripts/plot_electricity_production.py"
