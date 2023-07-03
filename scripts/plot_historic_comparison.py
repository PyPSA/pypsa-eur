#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Created on Mon Jul  3 12:50:26 2023.

@author: fabian
"""

import pandas as pd
import pypsa
from pypsa.statistics import get_bus_and_carrier

carrier_groups = {
    "Offshore Wind (AC)": "Offshore Wind",
    "Offshore Wind (DC)": "Offshore Wind",
    "Open-Cycle Gas": "Gas",
    "Combined-Cycle Gas": "Gas",
    "Reservoir & Dam": "Hydro",
    "Pumped Hydro Storage": "Hydro",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_statistics",
            simpl="",
            opts="Co2L-3h",
            clusters="37c",
            ll="v1.0",
        )

n = pypsa.Network(snakemake.input.network)
historic = pd.read_csv(
    snakemake.input.historic_electricity_generation, index_col=0, header=[0, 1]
)


simulated = n.statistics.dispatch(groupby=get_bus_and_carrier, aggregate_time=False).T
simulated = simulated[["Generator", "StorageUnit"]].droplevel(0, axis=1)
simulated = simulated.rename(columns=n.buses.country, level=0)
simulated = simulated.rename(carrier_groups, level=1)
simulated = simulated.groupby(axis=1, level=[0, 1]).sum()
