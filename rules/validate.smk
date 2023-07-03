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
        "../scripts/retrieve_electricity_production.py"


rule plot_electricity_production:
    input:
        network=RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        electricity_production="data/historical_electricity_production.csv",
    output:
        electricity_producion=RESULTS
        + "figures/validate_electricity_production_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.pdf",
    script:
        "scripts/plot_electricity_production.py"
