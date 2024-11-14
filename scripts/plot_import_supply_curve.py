# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates stacked bar charts of import shares per carrier.
"""

import logging

logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt
import pandas as pd
from _helpers import configure_logging

from plot_import_world_map import rename

TECH_MAPPING = {
    "hvdc-to-elec": "electricity (HVDC)",
    "pipeline-h2": "hydrogen (pipeline)",
    "shipping-lh2": "hydrogen (ship)",
    "shipping-lnh3": "ammonia (ship)",
    "shipping-lch4": "methane (ship)",
    "shipping-meoh": "methanol (ship)",
    "shipping-ftfuel": "Fischer-Tropsch (ship)",
    "shipping-hbi": "HBI (ship)",
    "shipping-steel": "steel (ship)",
}

import country_converter as coco
cc = coco.CountryConverter()

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_supply_curve",
            configfiles="config/config.20240826-z1.yaml",
        )

    configure_logging(snakemake)

    plt.style.use(["bmh", snakemake.input.rc])

    tech_colors = snakemake.config["plotting"]["tech_colors"]

    COLORS = {
        "wind": tech_colors["onwind"],
        "solar": tech_colors["solar"],
        "battery": tech_colors["battery"],
        "electrolysis": tech_colors["H2 Electrolysis"],
        "fuel storage": '#ffd4dc',
        "hydrogen conversion": tech_colors["Fischer-Tropsch"],
        "direct air capture": tech_colors["DAC"],
        "iron ore": "#4e4f55",
        "direct iron reduction": tech_colors["steel"],
        "electric arc furnace": "#8795a8",
        "evaporation/liquefaction": "#8487e8",
        "transport": "#e0ae75",
    }

    for esc, esc_nice_name in TECH_MAPPING.items():

        print(esc_nice_name)

        production = 100e6 if esc in ["shipping-steel", "shipping-hbi"] else 500e6

        df = pd.read_parquet(snakemake.input.imports).reset_index().query("scenario == 'default' and year == 2040 and category == 'cost' and esc == @esc")

        df.drop(["scenario", "year", "wacc", "esc", "category"], axis=1, inplace=True)

        df["subcategory"] = df["subcategory"].apply(rename)

        df = df.groupby(["exporter", "importer", "subcategory"]).value.sum().div(production).unstack("importer").min(axis=1).unstack("subcategory")

        cols = pd.Index(COLORS.keys()).intersection(df.columns)
        df = df.loc[df.sum(axis=1).sort_values(ascending=False).index, cols].reindex(COLORS.keys(), axis=1)

        split_index = df.index.str.split("-")
        suffix = split_index.str.get(1).fillna("").map(lambda x: f" ({x})" if x else x)
        prefix = pd.Index(cc.convert(split_index.str[0], to="short_name")).map(lambda x: "USA" if x == "United States" else x)
        df.index = prefix + suffix

        fig, ax = plt.subplots(figsize=(5, 11))

        df.plot.barh(ax=ax, stacked=True, color=COLORS, width=0.85)

        for i, (idx, row) in enumerate(df.iterrows()):
            sum_value = row.sum()
            ax.text(sum_value + 5, i, f"{sum_value:.1f}", va='center', ha='left', fontsize=8)

        handles, labels = ax.get_legend_handles_labels()
        handles.reverse()
        labels.reverse()
        ax.legend(handles, labels, title="", ncol=2, loc=(-0.55, 1.05))
        limit = 800 if esc in ["shipping-steel", "shipping-hbi"] else 200
        ax.set_xlim(0, limit + int(limit / 8))
        ax.set_xticks(range(0, limit + 1, int(limit / 4)))
        ax.set_xticks(range(int(limit / 8), limit + 1, int(limit / 8)), minor=True)
        ax.set_xlabel(esc_nice_name + (" [€/t]" if esc in ["shipping-steel", "shipping-hbi"] else " [€/MWh]"), fontsize=10)
        ax.set_ylabel("")
        ax.grid(False, axis="both")

        # Mirror the xticks on the top side of the axis
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ax.get_xticks())
        ax2.set_xticks(ax.get_xticks(minor=True), minor=True)
        ax2.set_xlabel(ax.get_xlabel(), fontsize=10)
        ax2.grid(False, axis="both")

        plt.tight_layout()

        for fn in snakemake.output[esc.replace("-", "_")]:
            plt.savefig(fn, bbox_inches="tight")
