# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Build afforestation CO2 sequestration potentials per network node.

Uses rotation-averaged NUTS-2 CO2 sequestration rates (tCO2/ha/y) from forest
yield tables, with monthly weights from FluxCom gross primary production, and
the CORINE land available per node. Network regions are intersected with the
NUTS-2 regions and the rates, carbon densities, rotation ages and monthly
weights are area-weighted per node.

Outputs: per node, the CO2 sequestration rate, potential [tCO2/y], carbon
density at rotation [tCO2/ha] and rotation age [y]; the node-level monthly
weights; the hourly seasonal profile; a map of the monthly rates.
"""

import logging

import geopandas
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas

from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)


MONTH_NAMES = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


def plot_afforestation_map(network, data_frame, monthly_weights_rows, output_png):
    """
    Save a 12-panel map of monthly afforestation CO2 sequestration rates
    [tCO2/(ha·month)] per clustered node.
    """
    gdf = network[["name", "geometry"]].set_index("name").to_crs(4326)
    gdf = gdf.join(data_frame)

    if monthly_weights_rows:
        mw = pandas.DataFrame(monthly_weights_rows).set_index("node")
        # monthly rate [tCO2/(ha·month)] = annual_rate × monthly_weight
        monthly_rates = mw[MONTH_NAMES].multiply(
            gdf["CO2 seq rate tCO2/(ha y)"], axis=0
        )
        gdf = gdf.join(monthly_rates)

        vmax = monthly_rates.values.max()
        fig, axes = plt.subplots(3, 4, figsize=(18, 10))
        axes = axes.ravel()
        cmap = "YlGn"

        for m_name, ax in zip(MONTH_NAMES, axes):
            gdf.plot(
                column=m_name,
                ax=ax,
                cmap=cmap,
                vmin=0,
                vmax=vmax,
                linewidth=0.3,
                edgecolor="white",
                missing_kwds={"color": "lightgrey"},
            )
            ax.set_title(m_name, fontsize=10)
            ax.set_xlim(-25, 45)
            ax.set_ylim(34, 72)
            ax.axis("off")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=0, vmax=vmax))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
        cbar.set_label("CO₂ sequestration rate (tCO₂ ha⁻¹ month⁻¹)", fontsize=9)
        fig.suptitle(
            "Monthly afforestation CO₂ sequestration rates per node\n"
            "(rotation-averaged MAI from Pilli et al. 2024, seasonalised with FluxCom GPP)",
            fontsize=11,
            y=1.01,
        )
    else:
        logger.warning("No monthly weights, skipping the afforestation map.")
        return

    fig.tight_layout()
    logger.info(f"Saving afforestation map to '{output_png}'")
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def build_afforestation_potentials(
    network_geojson,
    nuts2_geojson,
    corine_potentials_csv,
    afforestation_nuts_file,
    output_csv,
    monthly_weights_file=None,
    output_monthly_weights_csv=None,
    snapshot_config=None,
    output_seasonal_profile_csv=None,
    output_png=None,
):
    network = geopandas.read_file(network_geojson)
    corine_potentials = pandas.read_csv(corine_potentials_csv).set_index("node")

    logger.info("Calculating afforestation potentials based on growth rates.")
    nuts2 = geopandas.read_file(nuts2_geojson)
    nuts2_growth_rates = pandas.read_csv(afforestation_nuts_file).set_index("NUTS2")
    nuts2_monthly_weights = pandas.read_csv(monthly_weights_file, index_col="NUTS_ID")

    # harmonize NUTS2 indexes (EL → GR, UK → GB)
    def harmonize_idx(idx):
        return idx.map(
            lambda x: (
                ("GR" if x[:2] == "EL" else ("GB" if x[:2] == "UK" else x[:2])) + x[2:]
            )
        )

    nuts2["NUTS_ID"] = harmonize_idx(nuts2["NUTS_ID"])
    nuts2_growth_rates.index = harmonize_idx(nuts2_growth_rates.index)
    nuts2_monthly_weights.index = harmonize_idx(nuts2_monthly_weights.index)

    # Supplement nuts2 GeoDataFrame with pseudo-NUTS2 codes (e.g. RS00, AL00,
    # BA00, XK00) present in the rates CSV but absent from the Eurostat GeoJSON.
    # Use dissolved network-node geometry for the matching country prefix so
    # the subsequent overlay assigns the correct rate to every Balkan node.
    missing_codes = sorted(set(nuts2_growth_rates.index) - set(nuts2["NUTS_ID"]))
    if missing_codes:
        logger.info(
            f"Supplementing NUTS2 GDF with {len(missing_codes)} pseudo-regions "
            f"absent from GeoJSON: {missing_codes}"
        )
        network_4326 = network.to_crs(4326)
        extra_rows = []
        for code in missing_codes:
            country_prefix = code[:2]
            country_nodes = network_4326[
                network_4326["name"].str.startswith(country_prefix)
            ]
            if not country_nodes.empty:
                geom = country_nodes.geometry.union_all()
                extra_rows.append({"NUTS_ID": code, "geometry": geom})
                logger.info(
                    f"  {code}: geometry from {len(country_nodes)} node(s): "
                    f"{country_nodes['name'].tolist()}"
                )
            else:
                logger.warning(
                    f"  {code}: no network nodes found for country prefix '{country_prefix}' — skipping"
                )
        if extra_rows:
            extra_gdf = geopandas.GeoDataFrame(extra_rows, crs="EPSG:4326")
            nuts2 = geopandas.GeoDataFrame(
                pandas.concat([nuts2, extra_gdf], ignore_index=True), crs=nuts2.crs
            )

    # Build NUTS2 GeoDataFrame with rates and monthly weights joined.
    # Also join density and rotation_age for node-specific cost calculation.
    nuts2_gdf = nuts2[["NUTS_ID", "geometry"]].set_index("NUTS_ID")
    rate_cols = ["CO2 seq rate tCO2/(ha y)", "density tCO2/ha", "rotation_age_years"]
    nuts2_gdf = nuts2_gdf.join(nuts2_growth_rates[rate_cols])
    nuts2_gdf = nuts2_gdf.join(nuts2_monthly_weights[MONTH_NAMES])

    # Warn about NUTS2 regions missing monthly weights; fill with uniform 1/12 - it should not happen
    missing_w = nuts2_gdf[nuts2_gdf[MONTH_NAMES].isna().any(axis=1)].index.tolist()
    if missing_w:
        logger.warning(
            f"NUTS2 regions missing monthly weights (using uniform 1/12): {missing_w}"
        )
    nuts2_gdf[MONTH_NAMES] = nuts2_gdf[MONTH_NAMES].fillna(1.0 / 12)
    nuts2_gdf = nuts2_gdf.reset_index()  # NUTS_ID back as column for overlay

    # Project to EPSG:3035 for area calculations
    network_3035 = network.to_crs(3035)
    nuts2_3035 = nuts2_gdf.to_crs(3035)

    # Spatial overlay: one call replaces the manual double loop
    overlay = geopandas.overlay(
        network_3035[["name", "geometry"]],
        nuts2_3035,
        how="intersection",
        keep_geom_type=True,
    )

    # Share of node area covered by each NUTS2 region (weights for intensive quantities)
    node_areas = network_3035.set_index("name").geometry.area
    overlay["share"] = overlay.geometry.area / overlay["name"].map(node_areas)

    # CORINE available land per node [ha]
    overlay["total_area_ha"] = overlay["name"].map(
        corine_potentials["potential [sqkm]"] * 100
    )

    # Weighted contributions per intersection fragment
    overlay["potential_contrib"] = (
        overlay["total_area_ha"]
        * overlay["CO2 seq rate tCO2/(ha y)"]
        * overlay["share"]
    )
    overlay["density_contrib"] = overlay["density tCO2/ha"] * overlay["share"]
    overlay["rotage_contrib"] = overlay["rotation_age_years"] * overlay["share"]
    for m in MONTH_NAMES:
        overlay[f"w_{m}"] = overlay[m] * overlay["share"]

    # Aggregate per node
    w_cols = [f"w_{m}" for m in MONTH_NAMES]
    agg = overlay.groupby("name").agg(
        {
            "potential_contrib": "sum",
            "density_contrib": "sum",
            "rotage_contrib": "sum",
            "share": "sum",
            **{c: "sum" for c in w_cols},
        }
    )

    # Normalise monthly weights to sum to 1
    agg[w_cols] = agg[w_cols].div(agg[w_cols].sum(axis=1), axis=0)

    # Weighted-average growth rate [tCO2/(ha·y)]
    node_area_ha = corine_potentials["potential [sqkm]"] * 100
    agg["CO2 seq rate tCO2/(ha y)"] = agg["potential_contrib"] / node_area_ha.reindex(
        agg.index
    )

    # Area-weighted density and rotation age (normalised by actual share coverage)
    agg["density tCO2/ha"] = agg["density_contrib"] / agg["share"]
    agg["rotation_age_years"] = agg["rotage_contrib"] / agg["share"]

    for node_name, row in agg.iterrows():
        logger.info(
            f"Node '{node_name}': afforestation potential = "
            f"{row['potential_contrib']:.1f} tCO2/y (CO2 seq rate = "
            f"{row['CO2 seq rate tCO2/(ha y)']:.3f} tCO2/(ha y), density = "
            f"{row['density tCO2/ha']:.1f} tCO2/ha, T* = "
            f"{row['rotation_age_years']:.1f} yr)"
        )

    # Build data_frame with "node" as column (set_index called in shared save block below)
    agg.index.name = "node"
    out_cols = [
        "CO2 seq rate tCO2/(ha y)",
        "potential_contrib",
        "density tCO2/ha",
        "rotation_age_years",
    ]
    data_frame = (
        agg[out_cols]
        .rename(columns={"potential_contrib": "potential [tCO2/y]"})
        .reset_index()
    )

    # Node-level monthly weights for seasonal profile generation
    monthly_weights_rows = [
        {"node": node, **{m: agg.at[node, f"w_{m}"] for m in MONTH_NAMES}}
        for node in agg.index
    ]

    data_frame.set_index("node", inplace=True)
    logger.info(f"Saving afforestation potentials to '{output_csv}'")
    data_frame.to_csv(output_csv)

    mw_df = pandas.DataFrame(monthly_weights_rows).set_index("node")
    if output_monthly_weights_csv is not None:
        logger.info(
            f"Saving node-level monthly weights to '{output_monthly_weights_csv}'"
        )
        mw_df.to_csv(output_monthly_weights_csv)

    if output_seasonal_profile_csv is not None:
        if snapshot_config is not None:
            logger.info("Building hourly seasonal profile ...")
            snapshots = pandas.date_range(
                start=snapshot_config["start"],
                end=snapshot_config["end"],
                freq=snapshot_config.get("freq", "h"),
                inclusive="left",
            )
            T_year = len(snapshots)
            H_m = {m: int((snapshots.month == m).sum()) for m in range(1, 13)}
            profile = pandas.DataFrame(
                index=snapshots, columns=mw_df.index, dtype=float
            )
            for m_idx, m_name in enumerate(MONTH_NAMES, 1):
                mask = snapshots.month == m_idx
                h_m = H_m[m_idx]
                scale = T_year / h_m if h_m > 0 else 1.0 / 12
                profile.loc[mask] = (mw_df[m_name] * scale).values
            profile.index.name = "snapshot"
            logger.info(f"Saving seasonal profile to '{output_seasonal_profile_csv}'")
            profile.to_csv(output_seasonal_profile_csv)
        else:
            # no snapshot config: write empty placeholder
            pandas.DataFrame(columns=["snapshot"]).set_index("snapshot").to_csv(
                output_seasonal_profile_csv
            )

    if output_png is not None:
        plot_afforestation_map(network, data_frame, monthly_weights_rows, output_png)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_afforestation_potentials", clusters="39")

    configure_logging(snakemake)

    build_afforestation_potentials(
        network_geojson=snakemake.params["network_geojson"],
        nuts2_geojson=snakemake.params["nuts2_geojson"],
        corine_potentials_csv=snakemake.input[
            "afforestation_corine_potentials_csv_file"
        ],
        afforestation_nuts_file=snakemake.input["afforestation_nuts_file"],
        output_csv=snakemake.output["csv_file"],
        monthly_weights_file=snakemake.input["afforestation_monthly_weights_file"],
        output_monthly_weights_csv=snakemake.output["monthly_weights_csv_file"],
        snapshot_config=snakemake.params.get("snapshots"),
        output_seasonal_profile_csv=snakemake.output["seasonal_profile_csv_file"],
        output_png=snakemake.output["png_file"],
    )
