# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import pandas as pd
import pypsa
from _helpers import configure_logging, set_scenario_config
from entsoe import EntsoePandasClient
from entsoe.exceptions import InvalidBusinessParameterError, NoMatchingDataError
from requests import HTTPError

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_cross_border_flows")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    api_key = snakemake.config["private"]["keys"]["entsoe_api"]
    client = EntsoePandasClient(api_key=api_key)

    n = pypsa.Network(snakemake.input.network)
    start = pd.Timestamp(snakemake.params.snapshots["start"], tz="Europe/Brussels")
    end = pd.Timestamp(snakemake.params.snapshots["end"], tz="Europe/Brussels")

    branches = n.branches().query("carrier in ['AC', 'DC']")
    c = n.buses.country
    branch_countries = pd.concat([branches.bus0.map(c), branches.bus1.map(c)], axis=1)
    branch_countries = branch_countries.query("bus0 != bus1")
    branch_countries = branch_countries.apply(sorted, axis=1, result_type="broadcast")
    country_pairs = branch_countries.drop_duplicates().reset_index(drop=True)

    flows = []
    unavailable_borders = []
    for from_country, to_country in country_pairs.values:
        try:
            flow_directed = client.query_crossborder_flows(
                from_country, to_country, start=start, end=end
            )
            flow_reverse = client.query_crossborder_flows(
                to_country, from_country, start=start, end=end
            )
            flow = (flow_directed - flow_reverse).rename(
                f"{from_country} - {to_country}"
            )
            flow = flow.tz_localize(None).resample("1h").mean()
            flow = flow.loc[start.tz_localize(None) : end.tz_localize(None)]
            flows.append(flow)
        except (HTTPError, NoMatchingDataError, InvalidBusinessParameterError):
            unavailable_borders.append(f"{from_country}-{to_country}")

    if unavailable_borders:
        logger.warning(
            "Historical electricity cross-border flows for countries"
            f" {', '.join(unavailable_borders)} not available."
        )

    flows = pd.concat(flows, axis=1)
    flows.to_csv(snakemake.output[0])
