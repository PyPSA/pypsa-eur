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
