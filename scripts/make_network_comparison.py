# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Creates plots to compare transmission line lengths between two PyPSA networks.
Uses spatial intersection with country shapes for accurate country attribution.

"""

import logging

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
import pypsa
import seaborn as sns
from matplotlib.figure import Figure
from shapely import wkt
from shapely.geometry import LineString

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def thousands_formatter(x: float, pos: int) -> str:
    """Format axis values as thousands with 'k' suffix."""
    return f"{x * 1e-3:.0f}k"


def prepare_line_geometries(
    network: pypsa.Network, carrier: str = "AC"
) -> gpd.GeoDataFrame:
    """
    Prepare line geometries from PyPSA network.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network to extract lines from.
    carrier : str, default "AC"
        Carrier type - "AC" for lines, "DC" for links.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with line geometries and attributes.
    """
    if carrier == "AC":
        df = network.lines[["bus0", "bus1", "length", "num_parallel"]].copy()
        geom_source = network.lines
    else:  # DC
        df = network.links[network.links["carrier"] == "DC"].copy()
        if len(df) == 0:
            return gpd.GeoDataFrame(
                columns=["length", "num_parallel", "geometry"], crs="EPSG:4326"
            )
        df = df[["bus0", "bus1", "length"]].copy()
        df["num_parallel"] = 1
        geom_source = network.links[network.links["carrier"] == "DC"]

    # Use existing geometry if available
    if "geometry" in geom_source.columns and len(geom_source) > 0:
        geoms = geom_source["geometry"]
        # Handle WKT strings
        if isinstance(geoms.iloc[0], str):
            df["geometry"] = geoms.apply(wkt.loads)
        else:
            df["geometry"] = geoms.values
    else:
        # Fallback to straight line between buses
        logger.warning(
            f"No geometry column found for {carrier}, using straight-line approximation"
        )
        bus_coords = network.buses[["x", "y"]]
        df["geometry"] = df.apply(
            lambda row: LineString(
                [
                    (bus_coords.loc[row.bus0, "x"], bus_coords.loc[row.bus0, "y"]),
                    (bus_coords.loc[row.bus1, "x"], bus_coords.loc[row.bus1, "y"]),
                ]
            ),
            axis=1,
        )

    return gpd.GeoDataFrame(
        df[["length", "num_parallel"]], geometry=df["geometry"], crs="EPSG:4326"
    )


def get_lines_by_country(
    network: pypsa.Network, country_shapes: gpd.GeoDataFrame
) -> pd.DataFrame:
    """
    Calculate transmission line route and circuit lengths by country using spatial intersection.

    This function clips transmission lines at country boundaries to accurately attribute
    line lengths to countries, including cross-border lines. Uses spatial indexing and
    vectorized operations for performance.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network to analyze.
    country_shapes : gpd.GeoDataFrame
        GeoDataFrame with country geometries and 'name' column containing country codes.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['length_routes', 'length_circuits'] indexed by country.
    """

    # Prepare AC lines and DC links
    gdf_ac = prepare_line_geometries(network, carrier="AC")
    gdf_dc = prepare_line_geometries(network, carrier="DC")

    # Combine if both exist
    if len(gdf_ac) > 0 and len(gdf_dc) > 0:
        gdf_lines = pd.concat([gdf_ac, gdf_dc], ignore_index=True)
    elif len(gdf_ac) > 0:
        gdf_lines = gdf_ac
    elif len(gdf_dc) > 0:
        gdf_lines = gdf_dc
    else:
        logger.warning("No lines or links found in network")
        return pd.DataFrame(columns=["length_routes", "length_circuits"]).rename_axis(
            "country"
        )

    # Project to EPSG:3035 for accurate length calculations
    gdf_lines_proj = gdf_lines.to_crs("EPSG:3035")
    country_shapes_proj = country_shapes.to_crs("EPSG:3035")

    overlay = gpd.overlay(
        gdf_lines_proj,
        country_shapes_proj[["name", "geometry"]],
        how="intersection",
        keep_geom_type=True,  # Keep only LineString/MultiLineString
    )

    # Calculate lengths of clipped geometries
    overlay["clipped_length_km"] = overlay.geometry.length / 1000

    # Calculate route and circuit lengths
    overlay["length_routes"] = overlay["clipped_length_km"]
    overlay["length_circuits"] = overlay["clipped_length_km"] * overlay["num_parallel"]

    # Aggregate by country
    result = (
        overlay.groupby("name", observed=True)
        .agg({"length_routes": "sum", "length_circuits": "sum"})
        .rename_axis("country")
    )

    logger.info(f"Calculated line statistics for {len(result)} countries")
    logger.info(f"Total route length: {result['length_routes'].sum():.0f} km")
    logger.info(f"Total circuit length: {result['length_circuits'].sum():.0f} km")

    return result


def prepare_comparison_data(
    lines_incumbent: pd.DataFrame,
    lines_release: pd.DataFrame,
    countries: list,
    version: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Merge and prepare line length data for plotting.

    Parameters
    ----------
    lines_incumbent : pd.DataFrame
        Line lengths from baseline network, indexed by country.
    lines_release : pd.DataFrame
        Line lengths from new release network, indexed by country.
    countries : list
        Countries to include in comparison.
    version : str
        Version label for baseline network.

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame, pd.DataFrame)
        - routes_long: Long-format DataFrame for route lengths
        - circuits_long: Long-format DataFrame for circuit lengths
        - lines_merged: Wide-format merged DataFrame
    """
    # Merge datasets
    lines_merged = lines_incumbent.merge(
        lines_release, on="country", how="outer", suffixes=("_incumbent", "_release")
    )

    # Filter and sort by country
    lines_merged = (
        lines_merged.loc[lines_merged.index.intersection(countries)]
        .sort_index()
        .fillna(0)
    )

    # Label mapping
    label_map = {
        "incumbent": f"Incumbent network ({version})",
        "release": "New release",
    }

    def to_long_format(metric: str) -> pd.DataFrame:
        """Convert wide format to long format for a given metric."""
        col_prefix = f"length_{metric}_"
        return (
            lines_merged.reset_index()
            .melt(
                id_vars=["country"],
                value_vars=[f"{col_prefix}incumbent", f"{col_prefix}release"],
                var_name="parameter",
                value_name="length",
            )
            .assign(
                parameter=lambda df: df["parameter"]
                .str.replace(col_prefix, "", regex=False)
                .map(label_map)
            )
        )

    return to_long_format("routes"), to_long_format("circuits"), lines_merged


def plot_comparison(
    routes_data: pd.DataFrame, circuits_data: pd.DataFrame, fontsize: int = 10
) -> Figure:
    """
    Create comparison bar plots for route and circuit lengths.

    Parameters
    ----------
    routes_data : pd.DataFrame
        Long-format route length data with columns ['country', 'length', 'parameter'].
    circuits_data : pd.DataFrame
        Long-format circuit length data with columns ['country', 'length', 'parameter'].
    fontsize : int, default 10
        Font size for labels.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure.
    """
    palette = sns.color_palette("flare", n_colors=2)[::-1]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    # Create bar plots
    for ax, data, ylabel in zip(
        axes,
        [routes_data, circuits_data],
        ["a. Route length (km)", "b. Circuit length (km)"],
    ):
        sns.barplot(
            data=data, x="country", y="length", hue="parameter", ax=ax, palette=palette
        )
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(thousands_formatter))
        ax.legend().remove()

        # Add percentage change labels
        countries = data["country"].unique()
        for i, country in enumerate(countries):
            country_data = data[data["country"] == country]
            values = country_data.set_index("parameter")["length"]

            # Get comparison and release values
            incumbent_vals = values[values.index.str.contains("Incumbent")]
            release_vals = values[values.index.str.contains("New release")]

            if len(incumbent_vals) > 0 and len(release_vals) > 0:
                incumbent_val = incumbent_vals.iloc[0]
                release_val = release_vals.iloc[0]

                # Calculate percentage change
                if incumbent_val > 0:
                    pct_change = ((release_val - incumbent_val) / incumbent_val) * 100
                    y_pos = max(incumbent_val, release_val)

                    # Add text annotation
                    ax.text(
                        i,
                        y_pos,
                        f"{pct_change:+.1f}%",
                        ha="center",
                        va="bottom",
                        fontsize=fontsize - 2,
                        rotation=90,
                    )

    # Format x-axis
    axes[1].set_xlabel("Country", fontsize=fontsize)
    axes[1].tick_params(axis="x", rotation=45)

    # Adjust y-limits
    for ax in axes:
        ymax = ax.get_ylim()[1]
        ax.set_ylim(0, ymax * 1.15)

    # Add legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Route and circuit lengths",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=2,
        fontsize=fontsize,
        frameon=False,
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)

    return fig


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("make_network_comparison")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load networks
    logger.info("Loading networks...")
    n_release = pypsa.Network(snakemake.input.n_release)
    n_incumbent = pypsa.Network(snakemake.input.n_incumbent)
    countries = snakemake.params.countries
    version = snakemake.params.compare_to_version

    # Load and prepare country shapes
    logger.info(f"Loading country shapes from {snakemake.input.country_shapes}")
    country_shapes = gpd.read_file(snakemake.input.country_shapes)

    logger.info(f"Loading offshore regions from {snakemake.input.regions_offshore}")
    regions_offshore = gpd.read_file(snakemake.input.regions_offshore)[
        ["country", "geometry"]
    ]

    # Prepare offshore regions
    regions_offshore = regions_offshore.dissolve(by="country", as_index=False).rename(
        columns={"country": "name"}
    )
    regions_offshore = gpd.GeoDataFrame(regions_offshore)

    # Ensure country_shapes has the right columns
    country_shapes = gpd.GeoDataFrame(
        country_shapes[["name"]],
        geometry=country_shapes.geometry,
        crs=country_shapes.crs,
    )

    # Combine onshore and offshore regions
    regions = gpd.GeoDataFrame(
        pd.concat([country_shapes, regions_offshore], ignore_index=True),
        crs=country_shapes.crs,
    ).dissolve(by="name", as_index=False)

    # Filter to countries of interest
    if countries:
        regions = regions[regions["name"].isin(countries)]
        logger.info(f"Filtered to {len(regions)} countries/regions")

    # Calculate line statistics using spatial intersection
    logger.info("Calculating line statistics for incumbent network...")
    lines_incumbent = get_lines_by_country(n_incumbent, regions)

    logger.info("Calculating line statistics for release network...")
    lines_release = get_lines_by_country(n_release, regions)

    # Prepare comparison data
    routes_long, circuits_long, lines_merged = prepare_comparison_data(
        lines_incumbent,
        lines_release,
        countries,
        version,
    )

    # Calculate and log correlations
    pcc_routes = lines_merged["length_routes_incumbent"].corr(
        lines_merged["length_routes_release"]
    )
    pcc_circuits = lines_merged["length_circuits_incumbent"].corr(
        lines_merged["length_circuits_release"]
    )

    logger.info(f"Pearson correlation (route length): {pcc_routes:.4f}")
    logger.info(f"Pearson correlation (circuit length): {pcc_circuits:.4f}")

    # Create and save plot
    logger.info("Creating comparison plot...")
    fig = plot_comparison(routes_long, circuits_long, fontsize=10)

    logger.info(f"Exporting figure to {snakemake.output.lengths}")
    fig.savefig(snakemake.output.lengths, format="pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)
