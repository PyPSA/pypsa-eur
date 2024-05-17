# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Retrieve OSM data for the specified country using the overpass API and save it 
to the specified output files. Note that overpass requests are based on a fair 
use policy. `retrieve_osm_data` is meant to be used in a way that respects this 
policy by fetching the needed data once, only. 
"""

import json
import logging
import os
import requests
import time

from _helpers import configure_logging
logger = logging.getLogger(__name__)


# Function currently not needed - Kept for backup purposes to retrieve the OSM 
# area code if needed in the future
def _get_overpass_areas(countries):
    """
    Retrieve the OSM area codes for the specified country codes.
    
    Parameters
    ----------
    countries : str or list
        A single country code or a list of country codes for which the OSM area 
        codes should be retrieved.

    Returns
    -------
    dict
        A dictionary mapping country codes to their corresponding OSM area 
        codes.
    """

    # If a single country code is provided, convert it to a list
    if not isinstance(countries, list):
        countries = [countries]

    # Overpass API endpoint URL
    overpass_url = "https://overpass-api.de/api/interpreter"

    osm_areas = []
    for c in countries:
        # Overpass query to fetch the relation for the specified country code
        overpass_query = f"""
            [out:json];
            area["ISO3166-1"="{c}"];
            out;
        """

        # Send the request to Overpass API
        response = requests.post(overpass_url, data=overpass_query)

        try:
            # Parse the response
            data = response.json()

            # Check if the response contains any results
            if "elements" in data and len(data["elements"]) > 0:
                # Extract the area ID from the relation
                if c == "FR": # take second one for France
                    osm_area_id = data["elements"][1]["id"]
                else:
                    osm_area_id = data["elements"][0]["id"]
                osm_areas.append(f"area({osm_area_id})")
            else:
                # Print a warning if no results are found for the country code
                logger.info(f"No area code found for the specified country "
                            f"code: {c}. Omitted from the list.")
        except json.JSONDecodeError as e:
            logger.error(f"JSON decode error for country {c}: {e}")
            logger.debug(f"Response text: {response.text}")
    
    # Create a dictionary mapping country codes to their corresponding OSM area 
    # codes
    op_areas_dict = dict(zip(countries, osm_areas))
    
    return op_areas_dict
    

def retrieve_osm_data(
        country, 
        output,
        features=[
            "cables_way", 
            "lines_way", 
            "substations_way",
            "substations_relation",
            ]):
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
            "substations_way",
            "substations_relation",
            ].
    """
    # Overpass API endpoint URL
    overpass_url = "https://overpass-api.de/api/interpreter"

    # More features can in theory be retrieved that are currently not needed
    # to build a functioning network. The following power-related
    # features are supported:
    
    # features_dict= {
    #     'cables_way': 'way["power"="cable"]',
    #     'lines_way': 'way["power"="line"]',
    #     'substations_way': 'way["power"="substation"]',
    #     'substations_node': 'node["power"="substation"]',
    #     'transformers_way': 'way["power"="transformer"]',
    #     'transformers_node': 'node["power"="transformer"]',
    #     'route_relations': 'rel["route"="power"]["type"="route"]'
    # }

    features_dict= {
        'cables_way': 'way["power"="cable"]',
        'lines_way': 'way["power"="line"]',
        'substations_way': 'way["power"="substation"]',
        'substations_relation': 'relation["power"="substation"]',
    }

    wait_time = 5

    for f in features:
        if f not in features_dict:
            logger.info(f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}")
            raise ValueError(f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}")

        retries = 3
        for attempt in range(retries):
            logger.info(f" - Fetching OSM data for feature '{f}' in {country} (Attempt {attempt+1})...")

            # Build the overpass query
            op_area = f'area["ISO3166-1"="{country}"]'
            op_query = f'''
                [out:json];
                {op_area}->.searchArea;
                (
                {features_dict[f]}(area.searchArea);
                );
                out body geom;
            '''
            try:
                # Send the request
                response = requests.post(overpass_url, data = op_query)
                response.raise_for_status() # Raise HTTPError for bad responses
                data = response.json()

                filepath = output[f]
                parentfolder = os.path.dirname(filepath)
                if not os.path.exists(parentfolder):
                    os.makedirs(parentfolder)

                with open(filepath, mode = "w") as f:
                    json.dump(response.json(),f,indent=2)
                logger.info(" - Done.")
                break  # Exit the retry loop on success
            except (json.JSONDecodeError, requests.exceptions.RequestException) as e:
                logger.error(f"Error for feature '{f}' in country {country}: {e}")
                logger.debug(f"Response text: {response.text if response else 'No response'}")
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
                logger.error(f"Unexpected error for feature '{f}' in country {country}: {e}")
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

    # Retrieve the OSM data
    country = snakemake.wildcards.country
    output = snakemake.output

    retrieve_osm_data(country, output)