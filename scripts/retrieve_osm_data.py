# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve OSM data for the specified country using the overpass API and save it
to the specified output files.

Note that overpass requests are based on a fair
use policy. `retrieve_osm_data` is meant to be used in a way that respects this
policy by fetching the needed data once, only.
"""

import json
import logging
import os
import time

import requests
from _helpers import (  # set_scenario_config,; update_config_from_wildcards,; update_config_from_wildcards,
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def retrieve_osm_data(
    country,
    output,
    features=[
        "cables_way",
        "lines_way",
        "routes_relation",
        "substations_way",
        "substations_relation",
    ],
):
    """
    Retrieve OSM data for the specified country and save it to the specified
    output files.

    Parameters
    ----------
    country : str
        The country code for which the OSM data should be retrieved.
    output : dict
        A dictionary mapping feature names to the corresponding output file
        paths. Saving the OSM data to .json files.
    features : list, optional
        A list of OSM features to retrieve. The default is [
            "cables_way",
            "lines_way",
            "routes_relation",
            "substations_way",
            "substations_relation",
            ].
    """
    # Overpass API endpoint URL
    overpass_url = "https://overpass-api.de/api/interpreter"

    features_dict = {
        "cables_way": 'way["power"="cable"]',
        "lines_way": 'way["power"="line"]',
        "routes_relation": 'relation["route"="power"]',
        "substations_way": 'way["power"="substation"]',
        "substations_relation": 'relation["power"="substation"]',
    }

    wait_time = 5

    for f in features:
        if f not in features_dict:
            logger.info(
                f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}"
            )
            raise ValueError(
                f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}"
            )

        retries = 3
        for attempt in range(retries):
            logger.info(
                f" - Fetching OSM data for feature '{f}' in {country} (Attempt {attempt+1})..."
            )

            # Build the overpass query
            op_area = f'area["ISO3166-1"="{country}"]'
            op_query = f"""
                [out:json];
                {op_area}->.searchArea;
                (
                {features_dict[f]}(area.searchArea);
                );
                out body geom;
            """
            try:
                # Send the request
                response = requests.post(overpass_url, data=op_query)
                response.raise_for_status()  # Raise HTTPError for bad responses
                data = response.json()

                filepath = output[f]
                parentfolder = os.path.dirname(filepath)
                if not os.path.exists(parentfolder):
                    os.makedirs(parentfolder)

                with open(filepath, mode="w") as f:
                    json.dump(response.json(), f, indent=2)
                logger.info(" - Done.")
                break  # Exit the retry loop on success
            except (json.JSONDecodeError, requests.exceptions.RequestException) as e:
                logger.error(f"Error for feature '{f}' in country {country}: {e}")
                logger.debug(
                    f"Response text: {response.text if response else 'No response'}"
                )
                if attempt < retries - 1:
                    wait_time += 15
                    logger.info(f"Waiting {wait_time} seconds before retrying...")
                    time.sleep(wait_time)
                else:
                    logger.error(
                        f"Failed to retrieve data for feature '{f}' in country {country} after {retries} attempts."
                    )
            except Exception as e:
                # For now, catch any other exceptions and log them. Treat this
                # the same as a RequestException and try to run again two times.
                logger.error(
                    f"Unexpected error for feature '{f}' in country {country}: {e}"
                )
                if attempt < retries - 1:
                    wait_time += 10
                    logger.info(f"Waiting {wait_time} seconds before retrying...")
                    time.sleep(wait_time)
                else:
                    logger.error(
                        f"Failed to retrieve data for feature '{f}' in country {country} after {retries} attempts."
                    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_osm_data", country="BE")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Retrieve the OSM data
    country = snakemake.wildcards.country
    output = snakemake.output

    retrieve_osm_data(country, output)
