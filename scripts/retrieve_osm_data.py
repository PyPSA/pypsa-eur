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
import re
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
        "cables_way": [
            'way["power"="cable"]',
            'way["construction:power"="cable"]',
            'way["power"="construction"]["construction"="cable"]',
        ],
        "lines_way": [
            'way["power"="line"]',
            'way["construction:power"="line"]',
            'way["power"="construction"]["construction"="line"]',
        ],
        "routes_relation": [
            'relation["route"="power"]',
            'relation["power"="circuit"]',
            'relation["construction:power"="line"]',
            'relation["construction:power"="cable"]',
            'relation["power"="construction"]',
        ],
        "substations_way": [
            'way["power"="substation"]',
            'way["construction:power"="substation"]',
            'way["power"="construction"]["construction"="substation"]',
        ],
        "substations_relation": [
            'relation["power"="substation"]',
            'relation["construction:power"="substation"]',
            'relation["power"="construction"]["construction"="substation"]',
        ],
    }

    headers = {"User-Agent": user_agent}

    for f in features:
        if f not in features_dict:
            logger.info(
                f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}"
            )
            raise ValueError(
                f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}"
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

        wait_time = 30
        for attempt in range(max_tries):
            logger.info(
                f" - Fetching OSM data for feature '{f}' in {country} (Attempt {attempt + 1})..."
            )
            response = None
            try:
                # Client timeout slightly above the server-side query timeout
                response = requests.post(
                    url, data=op_query, headers=headers, timeout=timeout + 60
                )
                response.raise_for_status()  # Raise HTTPError for bad responses
                # Parse before opening the file, so that an error page or a
                # truncated response does not leave an empty output file behind
                data = response.json()
                # Overpass reports runtime errors during the query (e.g. timeout,
                # out of memory) as a remark next to incomplete results
                if "runtime error" in data.get("remark", ""):
                    raise ValueError(data["remark"])
            except (ValueError, requests.exceptions.RequestException) as e:
                # Overpass returns errors before the query as an HTML page with
                # status 200
                text = ""
                if response is not None and "html" in response.headers.get(
                    "content-type", ""
                ):
                    text = " ".join(re.sub(r"<[^>]+>", " ", response.text).split())
                logger.error(
                    f"Error for feature '{f}' in country {country}: {e}. "
                    f"Server response: {text[-300:]}"
                )
                if attempt == max_tries - 1:
                    raise RuntimeError(
                        f"Failed to retrieve data for feature '{f}' in country {country} after {max_tries} attempts."
                    ) from e
                # The same query is still running on the server, e.g. from an
                # aborted run, and is only released after its timeout
                if "duplicate_query" in text:
                    wait = timeout
                else:
                    wait = wait_time
                    wait_time *= 2
                logger.info(f"Waiting {wait} seconds before retrying...")
                time.sleep(wait)
                continue

            filepath = output[f]
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            with open(filepath, mode="w") as file:
                json.dump(data, file, indent=2)
            logger.info(" - Done.")
            break


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
