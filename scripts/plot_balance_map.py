# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create static energy balance maps for the defined carriers using`n.plot()`.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from packaging.version import Version, parse
from pypsa.plot import add_legend_lines, add_legend_patches, add_legend_semicircles
from pypsa.statistics import get_transmission_carriers

from scripts._helpers import (
    PYPSA_V1,
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import sanitize_carriers
from scripts.plot_power_network import load_projection

SEMICIRCLE_CORRECTION_FACTOR = 2 if parse(pypsa.__version__) <= Version("0.33.2") else 1

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_balance_map",
            clusters="50",
            opts="",
            sector_opts="",
            planning_horizons="2050",
            carrier="H2",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)
    sanitize_carriers(n, snakemake.config)
    pypsa.set_option("params.statistics.round", 8)
    pypsa.set_option("params.statistics.drop_zero", True)
    pypsa.set_option("params.statistics.nice_names", False)

    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    config = snakemake.params.plotting
    carrier = snakemake.wildcards.carrier
    settings = snakemake.params.settings
    carrier = carrier.replace(
        "_", " "
    )  # needed for slurm environment where [space] is not allowed

    # fill empty colors or "" with light grey
    mask = n.carriers.color.isna() | n.carriers.color.eq("")
    n.carriers["color"] = n.carriers.color.mask(mask, "lightgrey")

    # set EU location with location from config
    eu_location = config["eu_node_location"]
    n.buses.loc["EU", ["x", "y"]] = eu_location["x"], eu_location["y"]

    # get balance map plotting parameters
    boundaries = config["map"]["boundaries"]
    unit_conversion = settings["unit_conversion"]
    branch_color = settings.get("branch_color") or "darkseagreen"

    if carrier not in n.buses.carrier.unique():
        raise ValueError(
            f"Carrier {carrier} is not in the network. Remove from configuration `plotting: balance_map: bus_carriers`."
        )

    # for plotting change bus to location
    n.buses["location"] = n.buses["location"].replace("", "EU").fillna("EU")

    # set location of buses to EU if location is empty and set x and y coordinates to bus location
    n.buses["x"] = n.buses.location.map(n.buses.x)
    n.buses["y"] = n.buses.location.map(n.buses.y)

    # bus_size according to energy balance of bus carrier
    eb = n.statistics.energy_balance(bus_carrier=carrier, groupby=["bus", "carrier"])

    # remove energy balance of transmission carriers which relate to losses
    transmission_carriers = get_transmission_carriers(n, bus_carrier=carrier).rename(
        {"name": "carrier"}
    )
    components = transmission_carriers.unique("component")
    carriers = transmission_carriers.unique("carrier")

    # only carriers that are also in the energy balance
    carriers_in_eb = carriers[carriers.isin(eb.index.get_level_values("carrier"))]

    eb.loc[components] = eb.loc[components].drop(index=carriers_in_eb, level="carrier")
    eb = eb.dropna()
    bus_size = eb.groupby(level=["bus", "carrier"]).sum().div(unit_conversion)
    bus_size = bus_size.sort_values(ascending=False)

    # Get colors for carriers
    n.carriers.update({"color": snakemake.params.plotting["tech_colors"]})
    carrier_colors = n.carriers.color.copy().replace("", "grey")

    colors = (
        bus_size.index.get_level_values("carrier")
        .unique()
        .to_series()
        .map(carrier_colors)
    )

    # line and links widths according to optimal capacity
    flow = n.statistics.transmission(groupby=False, bus_carrier=carrier).div(
        unit_conversion
    )

    if not flow.empty:
        flow_reversed_mask = flow.index.get_level_values(1).str.contains("reversed")
        flow_reversed = flow[flow_reversed_mask].rename(
            lambda x: x.replace("-reversed", "")
        )
        flow = flow[~flow_reversed_mask].subtract(flow_reversed, fill_value=0)

    # if there are not lines or links for the bus carrier, use fallback for plotting
    fallback = pd.Series()
    line_width = flow.get("Line", fallback).abs()
    link_width = flow.get("Link", fallback).abs()

    # define maximal size of buses and branch width
    bus_size_factor = settings["bus_factor"]
    branch_width_factor = settings["branch_factor"]
    flow_size_factor = settings["flow_factor"]

    # get prices per region as colormap
    buses = n.buses.query("carrier in @carrier").index
    weights = n.snapshot_weightings.generators
    prices = weights @ n.buses_t.marginal_price[buses] / weights.sum()
    level = "name" if PYPSA_V1 else "Bus"
    price = prices.rename(n.buses.location).groupby(level=level).mean()

    if carrier == "co2 stored" and "CO2Limit" in n.global_constraints.index:
        co2_price = n.global_constraints.loc["CO2Limit", "mu"]
        price = price - co2_price

    # if only one price is available, use this price for all regions
    if price.size == 1:
        regions["price"] = price.values[0]
        shift = round(abs(price.values[0]) / 20, 0)
    else:
        regions["price"] = price.reindex(regions.index).fillna(0)
        shift = 0

    vmin, vmax = regions.price.min() - shift, regions.price.max() + shift
    if settings["vmin"] is not None:
        vmin = settings["vmin"]
    if settings["vmax"] is not None:
        vmax = settings["vmax"]

    crs = load_projection(snakemake.params.plotting)

    fig, ax = plt.subplots(
        figsize=(5, 6.5),
        subplot_kw={"projection": crs},
        layout="constrained",
    )

    line_flow = flow.get("Line")
    link_flow = flow.get("Link")
    transformer_flow = flow.get("Transformer")

    n.plot(
        bus_size=bus_size * bus_size_factor,
        bus_color=colors,
        bus_split_circle=True,
        line_width=line_width * branch_width_factor,
        link_width=link_width * branch_width_factor,
        line_flow=line_flow * flow_size_factor if line_flow is not None else None,
        link_flow=link_flow * flow_size_factor if link_flow is not None else None,
        link_color=branch_color,
        transformer_flow=transformer_flow * flow_size_factor
        if transformer_flow is not None
        else None,
        ax=ax,
        margin=0.2,
        geomap_color={"border": "darkgrey", "coastline": "darkgrey"},
        geomap=True,
        boundaries=boundaries,
    )

    regions.to_crs(crs.proj4_init).plot(
        ax=ax,
        column="price",
        cmap=settings["cmap"],
        vmin=vmin,
        vmax=vmax,
        edgecolor="None",
        linewidth=0,
    )

    ax.set_title(carrier)

    # Add colorbar
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=settings["cmap"], norm=norm)
    price_unit = settings["region_unit"]
    cbr = fig.colorbar(
        sm,
        ax=ax,
        label=f"Average Marginal Price [{price_unit}]",
        shrink=0.95,
        pad=0.03,
        aspect=50,
        orientation="horizontal",
    )
    cbr.outline.set_edgecolor("None")

    # add legend
    legend_kwargs = {
        "loc": "upper left",
        "frameon": False,
        "alignment": "left",
        "title_fontproperties": {"weight": "bold"},
    }

    pad = 0.18
    n.carriers.loc["", "color"] = "None"

    # Get lists for supply and consumption carriers
    pos_carriers = bus_size[bus_size > 0].index.unique("carrier")
    neg_carriers = bus_size[bus_size < 0].index.unique("carrier")

    # Determine larger total absolute value for supply and consumption for a carrier if carrier exists as both supply and consumption
    common_carriers = pos_carriers.intersection(neg_carriers)

    def get_total_abs(carrier, sign):
        values = bus_size.loc[:, carrier]
        return values[values * sign > 0].abs().sum()

    supp_carriers = sorted(
        set(pos_carriers) - set(common_carriers)
        | {c for c in common_carriers if get_total_abs(c, 1) >= get_total_abs(c, -1)}
    )
    cons_carriers = sorted(
        set(neg_carriers) - set(common_carriers)
        | {c for c in common_carriers if get_total_abs(c, 1) < get_total_abs(c, -1)}
    )

    # Add supply carriers
    add_legend_patches(
        ax,
        n.carriers.color[supp_carriers],
        supp_carriers,
        legend_kw={
            "bbox_to_anchor": (0, -pad),
            "ncol": 1,
            "title": "Supply",
            **legend_kwargs,
        },
    )

    # Add consumption carriers
    add_legend_patches(
        ax,
        n.carriers.color[cons_carriers],
        cons_carriers,
        legend_kw={
            "bbox_to_anchor": (0.5, -pad),
            "ncol": 1,
            "title": "Consumption",
            **legend_kwargs,
        },
    )

    # Add bus legend
    legend_bus_size = settings["bus_sizes"]
    carrier_unit = settings["unit"]
    if legend_bus_size is not None:
        add_legend_semicircles(
            ax,
            [
                s * bus_size_factor * SEMICIRCLE_CORRECTION_FACTOR
                for s in legend_bus_size
            ],
            [f"{s} {carrier_unit}" for s in legend_bus_size],
            patch_kw={"color": "#666"},
            legend_kw={
                "bbox_to_anchor": (0, 1),
                **legend_kwargs,
            },
        )

    # Add branch legend
    legend_branch_sizes = settings["branch_sizes"]
    if legend_branch_sizes is not None:
        add_legend_lines(
            ax,
            [s * branch_width_factor for s in legend_branch_sizes],
            [f"{s} {carrier_unit}" for s in legend_branch_sizes],
            patch_kw={"color": "#666"},
            legend_kw={"bbox_to_anchor": (0.25, 1), **legend_kwargs},
        )

    fig.savefig(
        snakemake.output[0],
        dpi=400,
        bbox_inches="tight",
    )
