# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
TODO To fill later
"""

import json
import logging
# import overpass as op
import os
import requests
import time

from _helpers import configure_logging
logger = logging.getLogger(__name__)


def _get_overpass_areas(countries):
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
            logger.info(f"No area code found for the specified country code: {c}. Ommitted from the list.")
    
    # Create a dictionary mapping country codes to their corresponding OSM area codes
    op_areas_dict = dict(zip(countries, osm_areas))
    
    return op_areas_dict
    

def retrieve_osm_data(
        country, 
        output,
        features=[
            "cables_way", 
            "lines_way", 
            "substations_way",
            "substations_node",
            "transformers_way",
            "transformers_node",
            "relations",
            ]):
    
    op_area = _get_overpass_areas(country)

    # Overpass API endpoint URL
    overpass_url = "https://overpass-api.de/api/interpreter"

    features_dict= {
        'cables_way': 'way["power"="cable"]',
        'lines_way': 'way["power"="line"]',
        'substations_way': 'way["power"="substation"]',
        'substations_node': 'node["power"="substation"]',
        'transformers_way': 'way["power"="transformer"]',
        'transformers_node': 'node["power"="transformer"]',
        'relations': 'rel["route"="power"]["type"="route"]'
    }

    for f in features:
        if f not in features_dict:
            raise ValueError(f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}")
            logger.info(f"Invalid feature: {f}. Supported features: {list(features_dict.keys())}")

        logger.info(f" - Fetching OSM data for feature '{f}' in {country}...")
        # Build the overpass query
        op_query = f'''
            [out:json];
            {op_area[country]}->.searchArea;
            (
            {features_dict[f]}(area.searchArea);
            );
            out body geom;
        '''

        # Send the request
        response = requests.post(overpass_url, data = op_query)
        # response = op.API(timeout=300).get(op_query) # returns data in geojson format. Timeout (max.) set to 300s

        filepath = output[f]
        parentfolder = os.path.dirname(filepath)
        if not os.path.exists(parentfolder):
            # Create the folder and its parent directories if they don't exist
            os.makedirs(parentfolder)

        with open(filepath, mode = "w") as f:
            # geojson.dump(response,f,indent=2)
            json.dump(response.json(),f,indent=2)
        logger.info(" - Done.")
        # time.sleep(5) 


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_osm_data", country="BE")
    
    configure_logging(snakemake)

    # Retrieve the OSM data
    country = snakemake.wildcards.country
    output = snakemake.output

    # Wait 5 seconds before fetching the OSM data to prevent too many requests error
    # TODO pypsa-eur: Add try catch to implement this only when needed
    logger.info(f"Waiting 5 seconds... Retrieving OSM data for {country}:")
    time.sleep(5) 
    retrieve_osm_data(country, output)