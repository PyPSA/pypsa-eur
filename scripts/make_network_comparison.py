# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Creates plots to compare transmission line lengths between two PyPSA networks.
"""

import logging

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
import pypsa
import seaborn as sns

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def thousands_formatter(x, pos):
    """Format axis values as thousands with 'k' suffix."""
    return f"{x * 1e-3:.0f}k"


def get_lines_by_country(network):
    """
    Calculate domestic transmission line route and circuit lengths by country.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network to analyze.

    Returns
    -------
    pd.DataFrame
        DataFrame with route and circuit lengths indexed by country.
    """
    # Get bus countries for both endpoints — lines
    lines = (
        network.lines[["bus0", "bus1", "length", "num_parallel"]]
        .merge(network.buses[["country"]], left_on="bus0", right_index=True)
        .merge(
            network.buses[["country"]],
            left_on="bus1",
            right_index=True,
            suffixes=("0", "1"),
        )
    )

    # Get bus countries for both endpoints — links
    links = (
        network.links.loc[network.links["carrier"] == "DC", ["bus0", "bus1", "length"]]
        .merge(network.buses[["country"]], left_on="bus0", right_index=True)
        .merge(
            network.buses[["country"]],
            left_on="bus1",
            right_index=True,
            suffixes=("0", "1"),
        )
    )
    links["num_parallel"] = 1  # Each link counts as one circuit

    # Combine lines and links
    combined = pd.concat([lines, links], ignore_index=True)

    # Aggregate parallel lines/links between same bus pairs
    combined = (
        combined.groupby(["bus0", "bus1", "country0", "country1"], observed=True)
        .agg({"length": "max", "num_parallel": "sum"})
        .reset_index()
    )

    # Keep only domestic (both endpoints in same country)
    combined = combined[combined["country0"] == combined["country1"]].copy()

    # Calculate route and circuit lengths
    combined["length_routes"] = combined["length"]
    combined["length_circuits"] = combined["num_parallel"] * combined["length"]

    # Sum by country
    return (
        combined.groupby("country0", observed=True)
        .agg({"length_routes": "sum", "length_circuits": "sum"})
        .rename_axis("country")
    )


def prepare_comparison_data(lines_incumbent, lines_release, countries, version):
    """
    Merge and prepare line length data for plotting.

    Parameters
    ----------
    lines_incumbent : pd.DataFrame
        Line lengths from baseline network
    lines_release : pd.DataFrame
        Line lengths from new release network
    countries : list
        Countries to include in comparison
    version : str
        Version label for baseline network

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame, pd.DataFrame)
        Long-format DataFrames for routes and circuits, and merged data
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

    # Single label mapping for all columns
    label_map = {
        "incumbent": f"Incumbent network ({version})",
        "release": "New release",
    }

    def to_long_format(metric):
        """Convert wide format to long format for a given metric (routes or circuits)."""
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
                .str.replace(col_prefix, "")
                .map(label_map)
            )
        )

    return to_long_format("routes"), to_long_format("circuits"), lines_merged


def plot_comparison(routes_data, circuits_data, fontsize=10):
    """
    Create comparison bar plots for route and circuit lengths.

    Parameters
    ----------
    routes_data : pd.DataFrame
        Long-format route length data
    circuits_data : pd.DataFrame
        Long-format circuit length data
    fontsize : int
        Font size for labels

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
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

            # Get comparison and release values (order depends on label mapping)
            comparison_label = [k for k in values.index if "Incumbent network" in k]
            release_label = [k for k in values.index if "New release" in k]

            if comparison_label and release_label:
                comparison_val = values[comparison_label[0]]
                release_val = values[release_label[0]]

                # Calculate percentage change
                if comparison_val > 0:
                    pct_change = ((release_val - comparison_val) / comparison_val) * 100

                    # Position label at top of the taller bar
                    y_pos = max(comparison_val, release_val)

                    # Format with sign
                    label = f"{pct_change:+.1f}%"

                    # Add text annotation
                    ax.text(
                        i,
                        y_pos,
                        label,
                        ha="center",
                        va="bottom",
                        fontsize=fontsize - 2,
                        rotation=90,
                    )

    # Format x-axis only on bottom plot
    axes[1].set_xlabel("Country", fontsize=fontsize)
    axes[1].tick_params(axis="x", rotation=45)

    ymin = 0
    for ax in axes:
        ymax = ax.get_ylim()[1]
        ax.set_ylim(ymin, ymax * 1.15)

    # Add single legend below figure
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

        snakemake = mock_snakemake(
            "make_network_comparison",
            configfiles="config/config.osm-release.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load networks
    n_release = pypsa.Network(snakemake.input.n_release)
    n_incumbent = pypsa.Network(snakemake.input.n_incumbent)
    countries = snakemake.params.countries
    version = snakemake.params.compare_to_version

    # Calculate line statistics
    lines_incumbent = get_lines_by_country(n_incumbent)
    lines_release = get_lines_by_country(n_release)

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
    fig = plot_comparison(routes_long, circuits_long, fontsize=10)

    logger.info(f"Exporting figure to {snakemake.output.lengths}")
    fig.savefig(snakemake.output.lengths, format="pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)
