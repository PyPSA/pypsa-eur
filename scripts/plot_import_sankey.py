# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates Sankey charts of import flows.
"""

import logging

import country_converter as coco
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pypsa
from _helpers import configure_logging
from plotly.io import write_image
from plotly.subplots import make_subplots
import plotly.io as pio   

pio.kaleido.scope.mathjax = None

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

TECH_COLORS = {
    "import pipeline-h2": "#db8ccd",
    "import shipping-lh2": "#e0c5dc",
    "import shipping-lnh3": "#e2ed74",
    "import shipping-lch4": "#f7a572",
    "import shipping-ftfuel": "#93eda2",
    "import shipping-meoh": "#87d0e6",
    "import shipping-hbi": "#d3dceb",
    "import shipping-steel": "#94a4be",
    "import hvdc-to-elec": "#91856a",
}

NICE_NAMES = {
    "import hvdc-to-elec": "electricity",
    "import pipeline-h2": "hydrogen (pipeline)",
    "import shipping-lh2": "hydrogen (ship)",
    "import shipping-ftfuel": "Fischer-Tropsch",
    "import shipping-meoh": "methanol",
    "import shipping-lch4": "methane (ship)",
    "import shipping-lnh3": "ammonia",
    "import shipping-steel": "steel",
    "import shipping-hbi": "HBI",
}


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake


        snakemake = mock_snakemake(
            "plot_import_sankey",
            opts="",
            clusters="115",
            ll="vopt",
            sector_opts="imp+AC+H20.9+CH40.9+NH30.9+FT0.9+MeOH0.9+HBI0.9+St0.9",
            planning_horizons="2050",
            configfiles="config/config.20240826-z1.yaml",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    df = (
        n.snapshot_weightings.generators @ n.links_t.p0.T.groupby(
            [n.links.carrier, n.links.bus0.str.replace(" export", ""), n.links.bus1.str[:2]]
        ).sum().T.filter(like="import").div(1e6)
    )

    df = df.loc[~df.index.get_level_values("carrier").str.contains("infrastructure")]
    df.name = "value"
    df = df.reset_index(level="carrier").reset_index()
    df.rename(columns=dict(bus0="source", bus1="target", carrier="label"), inplace=True)
    df = df.loc[df["value"] > 1]

    df["source"] = cc.convert(df.source.str.split("-").str[0], to="short_name")
    df["target"] = cc.convert(
        df.target.str.split("-").str[0], to="short_name", not_found="Anywhere in Europe"
    )
    country_aggregations = {
        "Tunisia": "Maghreb",
        "Algeria": "Maghreb",
        "Western Sahara": "Maghreb",
        "Morocco": "Maghreb",
        "Chile": "South America",
        "Argentina": "South America",
        "Libya": "Mashreq",
        "Egypt": "Mashreq",
    }
    df["source"] = df["source"].replace(country_aggregations)


    aggregate = df.groupby("target").value.sum() < 10
    aggregate = aggregate[aggregate].index.tolist()
    df["target"] = df.target.map(lambda x: "Other" if x in aggregate else x)

    import_volumes = df.groupby("target").value.sum().round(1)
    export_volumes = df.groupby("source").value.sum().round(1)

    df["source"] = df["source"].map(
        lambda x: f"{x} ({export_volumes[x]} TWh)"
    )
    df["target"] = df["target"].map(
        lambda x: f"{x} ({import_volumes[x]} TWh)"
    )

    labels = np.unique(df[["source", "target"]])
    nodes = pd.Series(range(len(labels)), index=labels)

    link_colors = df["label"].map(TECH_COLORS)

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    fig.add_trace(
        go.Sankey(
            arrangement="snap",
            valuesuffix=" TWh",
            valueformat=".1f",
            node=dict(
                pad=4,
                thickness=10,
                label=labels,
                color="#bbb",
                line=dict(color="black", width=0.7),
            ),
            link=dict(
                source=df.source.map(nodes),
                target=df.target.map(nodes),
                value=df.value,
                label=df.label,
                color=link_colors,
            ),
        )
    )

    for label, color in TECH_COLORS.items():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=color),
                name=NICE_NAMES[label],
                legendgroup=NICE_NAMES[label],
                showlegend=True,
            ),
            secondary_y=False,
        )

    axis_kwargs = dict(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_layout(
        height=500,
        width=400,
        margin=dict(l=20, r=20, t=20, b=20),
        font=dict(family="Helvetica, Arial", color="black"),
        paper_bgcolor="rgba(255,255,255,0.1)",
        plot_bgcolor="rgba(255,255,255,0.1)",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        xaxis=axis_kwargs,
        yaxis=axis_kwargs,
        xaxis2=axis_kwargs,
        yaxis2=axis_kwargs,
    )

    write_image(fig, snakemake.output[0])
