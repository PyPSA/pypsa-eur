# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot choropleth regional demands.
"""

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from _helpers import ensure_output_dir_exists
from pypsa.descriptors import get_switchable_as_dense as as_dense

TITLES = {
    "electricity": "Electricity Demand [TWh/a]",
    "H2": "Hydrogen Demand [TWh/a]",
    "heat": "Heat Demand [TWh/a]",
    "solid biomass": "Solid Biomass Demand [TWh/a]",
    "gas": "Methane Demand [TWh/a]",
    "oil": "Oil Demand [TWh/a]",
    "methanol": "Methanol Demand [TWh/a]",
    "ammonia": "Ammonia Demand [TWh/a]",
    "process emission": "Process Emissions [MtCO2/a]",
    # "steel": "Steel Demand [Mt/a]",
    # "HVC": "HVC Demand [Mt/a]",
}

CMAPS = {
    "electricity": "Blues",
    "H2": "RdPu",
    "heat": "Reds",
    "solid biomass": "Greens",
    "gas": "Oranges",
    "oil": "Greys",
    "methanol": "BuGn",
    "process emission": "Greys",
    "ammonia": "Purples",
    # "steel": "Blues",
    # "HVC": "Oranges",
}

REGEX = {
    "electricity": r"(electricity|EV)",
    "H2": r"(H2|fuel cell)",
    "heat": r"heat",
    "solid biomass": r"biomass",
    "oil": r"(oil|naphtha|kerosene)",
    "gas": r"gas",
    "methanol": r"methanol",
    "process emission": r"emission",
    "ammonia": r"ammonia",
    # "steel": r"(DRI|EAF|steel)",
    # "HVC": r"HVC",
}

MAPPING = {
    "H2 for industry": "H2",
    "agriculture electricity": "electricity",
    "agriculture heat": "heat",
    "agriculture machinery oil": "oil",
    "electricity": "electricity",
    "gas for industry": "gas",
    "industry electricity": "electricity",
    "industry methanol": "methanol",
    "kerosene for aviation": "oil",
    "land transport EV": "EV",
    "low-temperature heat for industry": "heat",
    "naphtha for industry": "oil",
    "urban central heat": "heat",
    "residential rural heat": "heat",
    "residential urban decentral heat": "heat",
    "services rural heat": "heat",
    "services urban decentral heat": "heat",
    "solid biomass for industry": "",
    "oil emissions": "emissions",
    "process emissions": "emissions",
    "agriculture machinery oil emissions": "emissions",
    "industry methanol emissions": "emissions",
}


def get_oil_demand(industrial_demand, nodal_energy_totals):
    oil = [
        "total international aviation",
        "total domestic aviation",
        "total agriculture machinery",
    ]
    return industrial_demand["naphtha"] + nodal_energy_totals[oil].sum(axis=1)


def get_biomass_demand(industrial_demand):
    return industrial_demand["solid biomass"]


def get_methanol_demand(
    industrial_demand, nodal_energy_totals, shipping_demand, config
):
    efficiency = (
        config["sector"]["shipping_oil_efficiency"]
        / config["sector"]["shipping_methanol_efficiency"]
    )
    return (
        industrial_demand["methanol"]
        + nodal_energy_totals["total domestic navigation"]
        + shipping_demand * efficiency
    )


def get_ammonia_demand(industrial_demand, shipping_demand, config):
    efficiency = (
        config["sector"]["shipping_oil_efficiency"]
        / config["sector"]["shipping_ammonia_efficiency"]
    )
    return industrial_demand["ammonia"] + shipping_demand * efficiency


def get_process_emissions(industrial_demand):
    cols = ["process emission", "process emission from feedstock"]
    return industrial_demand[cols].sum(axis=1)


def get_demand_by_region(
    n, industrial_demand, nodal_energy_totals, shipping_demand, config
):
    demand = as_dense(n, "Load", "p_set").div(1e6)  # TWh
    demand_grouped = demand.groupby(
        [n.loads.carrier, n.loads.bus.map(n.buses.location)], axis=1
    ).sum()
    df = (n.snapshot_weightings.generators @ demand_grouped).unstack(level=0)
    df.drop(["", "EU"], inplace=True)
    df.dropna(axis=1, how="all", inplace=True)
    oil_demand = pd.DataFrame(
        {"oil": get_oil_demand(industrial_demand, nodal_energy_totals)}
    )

    methanol_demand = pd.DataFrame(
        {
            "methanol": get_methanol_demand(
                industrial_demand, nodal_energy_totals, shipping_demand, config
            )
        }
    )

    biomass_demand = pd.DataFrame({"biomass": get_biomass_demand(industrial_demand)})

    process_emissions_demand = pd.DataFrame(
        {"process emissions": get_process_emissions(industrial_demand)}
    )

    ammonia_demand = pd.DataFrame(
        {"ammonia": get_ammonia_demand(industrial_demand, shipping_demand, config)}
    )

    df = pd.concat(
        [
            df,
            oil_demand,
            methanol_demand,
            biomass_demand,
            process_emissions_demand,
            ammonia_demand,
        ],
        axis=1,
    )

    return df


def plot_regional_demands(df, geodf, carrier, dir="."):
    series = df.filter(regex=REGEX[carrier]).sum(axis=1)

    proj = ccrs.EqualEarth()
    geodf = geodf.to_crs(proj.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": proj})

    geodf.plot(
        ax=ax,
        column=series,
        cmap=CMAPS[carrier],
        linewidths=0,
        legend=True,
        legend_kwds={"label": TITLES[carrier], "shrink": 0.7},
    )

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    ax.axis("off")

    carrier_fn = carrier.replace("-", "_").replace(" ", "_")
    fn = f"map-{carrier_fn}"
    plt.savefig(dir + "/" + fn + ".png", bbox_inches="tight")
    plt.savefig(dir + "/" + fn + ".pdf", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_choropleth_demand",
            simpl="",
            clusters=100,
            ll="vopt",
            opts="",
            sector_opts="Co2L0-73SN-T-H-B-I-S-A",
            planning_horizons=2050,
            configfiles=["../../config/config.20230825.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    ensure_output_dir_exists(snakemake)
    dir = snakemake.output[0]

    n = pypsa.Network(snakemake.input.network)

    industrial_demand = pd.read_csv(snakemake.input.industrial_demand, index_col=[0, 1])

    shipping_demand = pd.read_csv(
        snakemake.input.shipping_demand, index_col=0
    ).squeeze()

    nodal_energy_totals = pd.read_csv(snakemake.input.nodal_energy_totals, index_col=0)

    options = snakemake.config["sector"]
    endogenous_sectors = []
    if options["endogenous_steel"]:
        endogenous_sectors += ["DRI + Electric arc"]
    sectors_b = ~industrial_demand.index.get_level_values("sector").isin(
        endogenous_sectors
    )
    industrial_demand = industrial_demand.loc[sectors_b].groupby(level=0).sum()

    df = get_demand_by_region(
        n, industrial_demand, nodal_energy_totals, shipping_demand, snakemake.config
    )

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    for carrier in TITLES.keys():
        plot_regional_demands(df, regions, carrier, dir=dir)
