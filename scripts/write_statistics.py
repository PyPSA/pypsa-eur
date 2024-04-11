#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import seaborn as sns
from _helpers import configure_logging
from pypsa.statistics import get_carrier


# grouperfunctions = hier schreiben und dann in statistics.
def groupby_country_and_carrier(n, c, nice_names=False):
    df = n.df(c)
    bus = "bus1" if "bus" not in n.df(c) else "bus"
    country = df[bus].map(n.buses.location).map(n.buses.country).rename("country")
    carrier = get_carrier(n, c, nice_names)
    return [country, carrier]


def call_with_handle(func, **kwargs):
    try:
        ds = func(**kwargs)
    except Exception as e:
        print(f"An error occurred: {e}")
        ds = pd.Series()
        pass
    return ds


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "save_statistics_csv",
            simpl="",
            ll="v1.5",
            clusters="5",
            opts="",
            sector_opts="24h-T-H-B-I-A-dist1",
            planning_horizons="2040",
            country="all",
            carrier="electricity",
        )
    # configure_logging(snakemake)
    config = snakemake.config

    n = pypsa.Network(snakemake.input.network)

    kwargs = {"nice_names": False}

    wildcards = dict(snakemake.wildcards)
    carrier = wildcards.pop("carrier")
    country = wildcards.pop("country")

    if carrier == "all":
        pass
    elif carrier in config["plotting"].get("carrier_groups", []):
        bus_carrier = config["plotting"]["carrier_groups"][carrier]
        kwargs["bus_carrier"] = bus_carrier
    elif n.buses.carrier.str.contains(carrier).any():
        if carrier not in n.buses.carrier.unique():
            carrier = [
                bus_carrier
                for bus_carrier in n.buses.carrier.unique()
                if carrier in bus_carrier
            ]
        bus_carrier = carrier
        kwargs["bus_carrier"] = bus_carrier
    else:
        raise ValueError(
            f"Carrier {carrier} not found in network or carrier group in config."
        )

    if country == "all" or not country:
        kwargs["groupby"] = get_carrier
    else:
        kwargs["groupby"] = groupby_country_and_carrier

    for output in snakemake.output.keys():
        if "touch" in output:
            continue
        if output == "energy_balance":
            supply = call_with_handle(n.statistics.supply, **kwargs)
            withdrawal = call_with_handle(n.statistics.withdrawal, **kwargs)
            ds = (
                pd.concat([supply, withdrawal.mul(-1)])
                .groupby(supply.index.names)
                .sum()
            )
            ds.attrs = supply.attrs
            ds.attrs["name"] = "Energy Balance"
        elif output == "total_cost":
            opex = call_with_handle(n.statistics.opex, **kwargs)
            capex = call_with_handle(n.statistics.capex, **kwargs)
            ds = opex.add(capex, fill_value=0)
            ds.attrs["name"] = "Total Cost"
        else:
            func = eval(f"n.statistics.{output}")
            ds = call_with_handle(func, **kwargs)

        if ds.empty:
            print(
                f"Empty series for {output} with bus carrier {bus_carrier} and country {country}."
            )
            pd.Series().to_csv(snakemake.output[output])
            continue
        if country and country != "all":
            ds = ds.xs(country, level="country")
        pd.Series(ds.attrs).to_csv(snakemake.output[output], header=False)
        ds.dropna().to_csv(snakemake.output[output], mode="a")
    # touch file
    with open(snakemake.output.csv_touch, "a"):
        pass
