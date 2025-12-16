# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create interactive energy balance maps for the defined carriers using `n.explore()`.
"""

import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pydeck as pdk
import pypsa
from pypsa.plot.maps.interactive import PydeckPlotter
from pypsa.statistics import get_transmission_carriers

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.add_electricity import sanitize_carriers

VALID_MAP_STYLES = PydeckPlotter.VALID_MAP_STYLES


def scalar_to_rgba(
    value: float,
    *,
    norm: mcolors.Normalize,
    cmap: mcolors.Colormap,
    alpha: float = 1.0,
) -> list[int]:
    """
    Map a scalar float value to an RGBA color encoded as 8-bit integers.

    Parameters
    ----------
    value : float
        Scalar input to map through the normalization and colormap.
    norm : matplotlib.colors.Normalize
        Normalization defining vmin and vmax used for scaling.
    cmap : matplotlib.colors.Colormap
        Colormap used to convert normalized values to RGBA colors.
    alpha : float, optional (default = 1.0)
        Opacity in the range [0, 1]. Overrides the colormap's alpha.

    Returns
    -------
    List[int]
        A list ``[R, G, B, A]`` where each channel is an integer in the 0–255 range.
    """

    # Clamp to normalization bounds
    p = max(norm.vmin, min(norm.vmax, value))

    # Convert to RGBA floats (0–1)
    r, g, b, _ = cmap(norm(p))

    # Clamp and apply alpha
    a = alpha if 0.0 <= alpha <= 1.0 else 1.0

    # Convert to 8-bit integers
    return [
        int(r * 255),
        int(g * 255),
        int(b * 255),
        int(a * 255),
    ]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_balance_map_interactive",
            clusters=50,
            opts="",
            sector_opts="",
            planning_horizons="2050",
            carrier="H2",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Interactive map settings
    settings = snakemake.params.settings
    unit_conversion = settings["unit_conversion"]
    cmap = settings["cmap"]
    region_alpha = settings["region_alpha"]
    region_unit = settings["region_unit"]
    branch_color = settings["branch_color"]
    arrow_size_factor = settings["arrow_size_factor"]
    bus_size_max = settings["bus_size_max"]
    branch_width_max = settings["branch_width_max"]
    map_style = settings.get("map_style")
    map_style = VALID_MAP_STYLES.get(map_style, "road")
    tooltip = settings["tooltip"]

    # Import
    n = pypsa.Network(snakemake.input.network)
    sanitize_carriers(n, snakemake.config)
    pypsa.options.params.statistics.round = 8
    pypsa.options.params.statistics.drop_zero = True
    pypsa.options.params.statistics.nice_names = False

    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    carrier = snakemake.wildcards.carrier
    carrier = carrier.replace("_", " ")

    # Fill missing carrier colors
    missing_color = "#808080"
    b_missing = n.carriers.query("color == '' or color.isnull()").index
    n.carriers.loc[b_missing, "color"] = missing_color

    transmission_carriers = get_transmission_carriers(n, bus_carrier=carrier).rename(
        {"name": "carrier"}
    )
    components = transmission_carriers.unique("component")
    carriers = transmission_carriers.unique("carrier")

    ### Pie charts
    eb = n.statistics.energy_balance(
        bus_carrier=carrier,
        groupby=["bus", "carrier"],
    )

    # Only carriers that are also in the energy balance
    carriers_in_eb = carriers[carriers.isin(eb.index.get_level_values("carrier"))]

    eb.loc[components] = eb.loc[components].drop(index=carriers_in_eb, level="carrier")
    eb = eb.dropna()
    bus_size = eb.groupby(level=["bus", "carrier"]).sum()

    # line and links widths according to optimal capacity
    flow = n.statistics.transmission(groupby=False, bus_carrier=carrier)
    if not flow.empty:
        flow_reversed_mask = flow.index.get_level_values(1).str.contains("reversed")
        flow_reversed = flow[flow_reversed_mask].rename(
            lambda x: x.replace("-reversed", "")
        )
        flow = flow[~flow_reversed_mask].subtract(flow_reversed, fill_value=0)

    # only line first index
    line_flow = (
        flow.loc[flow.index.get_level_values(0).str.contains("Line")]
        .copy()
        .droplevel(0)
    )
    link_flow = (
        flow.loc[flow.index.get_level_values(0).str.contains("Link")]
        .copy()
        .droplevel(0)
    )

    branch_components = ["Link"]
    if carrier == "AC":
        branch_components = ["Line", "Link"]

    ### Prices
    buses = n.buses.query("carrier in @carrier").index
    demand = (
        n.statistics.energy_balance(
            bus_carrier=carrier, aggregate_time=False, groupby=["bus", "carrier"]
        )
        .clip(lower=0)
        .groupby("bus")
        .sum()
        .reindex(buses)
        .rename(n.buses.location)
        .T
    )

    weights = n.snapshot_weightings.generators
    price = (
        weights
        @ n.buses_t.marginal_price.reindex(buses, axis=1).rename(
            n.buses.location, axis=1
        )
        / weights.sum()
    )

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

    # Map colors
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap)

    regions["color"] = regions["price"].apply(
        scalar_to_rgba,
        norm=norm,
        cmap=cmap,
        alpha=region_alpha,
    )

    # Create tooltips
    regions["tooltip_html"] = (
        "<b>"
        + regions.index
        + "</b><br>"
        + "<b>Weighted price:</b> "
        + regions["price"].round(2).astype(str)
        + " "
        + region_unit
    )
    # regions["tooltip_html"] = regions["price"].round(2).astype(str)
    # Create layer
    regions_layer = pdk.Layer(
        "GeoJsonLayer",
        regions,
        stroked=True,
        filled=True,
        get_fill_color="color",
        get_line_color=[255, 255, 255, 255],
        line_width_min_pixels=1,
        pickable=True,
        auto_highlight=True,
    )

    map = n.explore(
        branch_components=branch_components,
        bus_size=bus_size.div(unit_conversion),
        bus_split_circle=True,
        line_width=line_flow.div(unit_conversion),
        line_flow=line_flow.div(unit_conversion),
        line_color="rosybrown",
        link_width=link_flow.div(unit_conversion),
        link_flow=link_flow.div(unit_conversion),
        link_color=branch_color,
        arrow_size_factor=arrow_size_factor,
        tooltip=tooltip,
        auto_scale=True,
        branch_width_max=branch_width_max,
        bus_size_max=bus_size_max,
        map_style=map_style,
    )

    map.layers.insert(0, regions_layer)

    map.to_html(snakemake.output[0], offline=True)
