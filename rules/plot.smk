# SPDX-FileCopyrightText: 2023- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

rule plot_power_network_unclustered:
    input:
        network=RESOURCES + "networks/elec.nc",
        rc="matplotlibrc",
    output:
        multiext(RESOURCES + "graphics/power-network-unclustered", ".png", ".pdf")
    script:
        "../scripts/plot_power_network_unclustered.py"


rule plot_gas_network_unclustered:
    input:
        gas_network=RESOURCES + "gas_network.csv",
        gem=HTTP.remote(
            "https://globalenergymonitor.org/wp-content/uploads/2023/07/Europe-Gas-Tracker-2023-03-v3.xlsx",
            keep_local=True,
        ),
        entry="data/gas_network/scigrid-gas/data/IGGIELGN_BorderPoints.geojson",
        storage="data/gas_network/scigrid-gas/data/IGGIELGN_Storages.geojson",
        rc="matplotlibrc",
    output:
        multiext(RESOURCES + "graphics/gas-network-unclustered", ".png", ".pdf")
    script:
        "../scripts/plot_gas_network_unclustered.py"
