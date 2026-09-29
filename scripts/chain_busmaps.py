# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build combined busmap from simplify_network and cluster_network busmaps.

Chains the two busmaps to create a single mapping from base network buses
to clustered network buses.

Inputs
------
- ``resources/busmap_simplify_network.csv``: Mapping from base to simplified buses
- ``resources/busmap_cluster_network.csv``: Mapping from simplified to clustered buses

Outputs
-------
- ``resources/busmap.csv``: Combined mapping from base to clustered buses
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, sanitize_busmap, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("chain_busmaps")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    simplify_busmap = pd.read_csv(
        snakemake.input.busmap_simplify_network, index_col=0
    ).squeeze()
    cluster_busmap = pd.read_csv(
        snakemake.input.busmap_cluster_network, index_col=0
    ).squeeze()

    busmap = sanitize_busmap(simplify_busmap.map(cluster_busmap))
    busmap.to_csv(snakemake.output[0])

    logger.info(f"Combined busmap with {len(busmap)} entries saved.")
