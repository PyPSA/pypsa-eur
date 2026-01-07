# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
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
    url="https://overpass-api.de/api/interpreter",
    max_tries=3,
    timeout=600,
    user_agent="",
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
    url : str, optional
        The URL of the overpass API endpoint. The default is
        "https://overpass-api.de/api/interpreter".
    max_tries : int, optional
        The maximum number of attempts to retrieve the data in case of failure. The
        default is 3.
    timeout : int, optional
        The timeout in seconds for the overpass API requests. The default is 600.
    user_agent : str
        The User-Agent string to include in the request headers for fair use policy compliance.

    """

    features_dict = {
        "cables_way": ['way["power"="cable"]'],
        "lines_way": ['way["power"="line"]'],
        "routes_relation": ['relation["route"="power"]', 'relation["power"="circuit"]'],
        "substations_way": ['way["power"="substation"]'],
        "substations_relation": ['relation["power"="substation"]'],
    }

    wait_time = 5

    headers = {"User-Agent": user_agent}

    for f in features:
        if f not in features_dict:
            logger.info(
                f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}"
            )
            raise ValueError(
                f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}"
            )

        max_tries = 3
        for attempt in range(max_tries):
            logger.info(
                f" - Fetching OSM data for feature '{f}' in {country} (Attempt {attempt + 1})..."
            )

            # Build the overpass query
            op_area = f'area["ISO3166-1"="{country}"]'
            op_query = f"""
                [out:json][timeout:{timeout}];
                {op_area}->.searchArea;
                (
                {" ".join(f"{i}(area.searchArea);" for i in features_dict[f])}
                );
                out body geom;
            """
            try:
                # Send the request
                response = requests.post(url, data=op_query, headers=headers)
                response.raise_for_status()  # Raise HTTPError for bad responses

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
                if attempt < max_tries - 1:
                    wait_time += 15
                    logger.info(f"Waiting {wait_time} seconds before retrying...")
                    time.sleep(wait_time)
                else:
                    logger.error(
                        f"Failed to retrieve data for feature '{f}' in country {country} after {max_tries} attempts."
                    )
            except Exception as e:
                # For now, catch any other exceptions and log them. Treat this
                # the same as a RequestException and try to run again two times.
                logger.error(
                    f"Unexpected error for feature '{f}' in country {country}: {e}"
                )
                if attempt < max_tries - 1:
                    wait_time += 10
                    logger.info(f"Waiting {wait_time} seconds before retrying...")
                    time.sleep(wait_time)
                else:
                    logger.error(
                        f"Failed to retrieve data for feature '{f}' in country {country} after {max_tries} attempts."
                    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "retrieve_osm_data_raw",
            country="BE",
        )

    overpass_api = snakemake.params.overpass_api
    url = overpass_api["url"]
    max_tries = overpass_api["max_tries"]
    timeout = overpass_api["timeout"]

    # Build User-Agent header
    ua_cfg = overpass_api["user_agent"]
    project = ua_cfg["project_name"]
    email = ua_cfg["email"]
    website = ua_cfg["website"]

    user_agent = f"{project} (Contact: {email}; Website: {website})"

    # Retrieve the OSM data
    country = snakemake.wildcards.country
    output = snakemake.output

    retrieve_osm_data(
        country,
        output,
        features=[
            "cables_way",
            "lines_way",
            "routes_relation",
            "substations_way",
            "substations_relation",
        ],
        url=url,
        max_tries=max_tries,
        timeout=timeout,
        user_agent=user_agent,
    )
