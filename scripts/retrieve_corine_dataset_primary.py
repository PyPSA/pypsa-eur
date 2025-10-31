# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
To download CORINE dataset from the primary data source - https://land.copernicus.eu/en/products/corine-land-cover/clc-2012

Usage Instructions:
    1. Login using EU login at https://land.copernicus.eu/user/login and create an API key
    2. Copy API key into the config.default.yaml -> (save from portal)
    #   secrets:
    #       corine: ''
"""

import json
import logging
import time
from json.decoder import JSONDecodeError
from pathlib import Path
from shutil import copy2, unpack_archive

import jwt
import requests

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_access_token(apikey):
    # Login using EU login at https://land.copernicus.eu/user/login and create an API key
    # Copy API key into the config.default.yaml -> (save from portal)
    #   secrets:
    #       corine: ''
    try:
        service_key = json.loads(apikey)
        private_key = service_key["private_key"].encode("utf-8")
        claim_set = {
            "iss": service_key["client_id"],
            "sub": service_key["user_id"],
            "aud": service_key["token_uri"],
            "iat": int(time.time()),
            "exp": int(time.time() + 3600),  # max 1 hour
        }
        assertion = jwt.encode(claim_set, private_key, algorithm="RS256")

        token_request = requests.post(
            service_key["token_uri"],
            headers={
                "Accept": "application/json",
                "Content-Type": "application/x-www-form-urlencoded",
            },
            data={
                "grant_type": "urn:ietf:params:oauth:grant-type:jwt-bearer",
                "assertion": assertion,
            },
        )
        token_request.raise_for_status()
        data = token_request.json()
        access_token = data.get("access_token")
    except JSONDecodeError as e:
        raise ValueError(
            "Missing or invalid access_token for corine. Check usage instructions in the scripts/retrieve_corine_dataset_primary.py script before proceeding"
        ) from e

    return access_token


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_corine_dataset_primary")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    apikey = snakemake.params["apikey"]
    output_zip_file = snakemake.output["zip"]
    tif_file = snakemake.output["tif_file"]
    access_token = load_access_token(apikey)

    if access_token:
        HEADERS = {
            "Accept": "application/json",
            "Content-Type": "application/json",
            "Authorization": f"Bearer {access_token}",
        }

        requests.post(
            "https://land.copernicus.eu/api/@datarequest_post",
            headers=HEADERS,
            json={
                "Datasets": [
                    {
                        "DatasetID": "a5ee71470be04d66bcff498f94ceb5dc",
                        "FileID": "9f992dad-d129-408e-be66-821adbd52a46",
                    }
                ]
            },
        )

        status_url = (
            "https://land.copernicus.eu/api/@datarequest_search?status=Finished_ok"
        )

        download_link = None

        logger.info("Waiting for the download to be prepared...")

        # Wait up to 20 minutes, checking every 60 seconds
        for attempt in range(20):
            time.sleep(60)  # avoid 429 Too Many Requests
            response = requests.get(status_url, headers=HEADERS)

            if response.status_code != 200:
                logger.info(
                    f"Attempt {attempt + 1}: Status check failed — {response.status_code}"
                )
                continue

            results = response.json()

            for result in results:
                if results[result].get("Status") == "Finished_ok":
                    download_link = results[result].get("DownloadURL")

                    if download_link:
                        file_response = requests.get(
                            download_link,
                            headers={
                                "Authorization": f"Bearer {access_token}",
                                "User-Agent": "curl/7.81.0",
                            },
                            allow_redirects=False,
                        )
                        if file_response.status_code == 200:
                            with open(output_zip_file, "wb") as f:
                                f.write(file_response.content)
                            logger.info(f"File saved as {output_zip_file}")

                            output_folder = Path(output_zip_file).parent
                            unpack_archive(output_zip_file, output_folder)

                            # unpack the actual dataset inside the downloaded zip -
                            # with new versions, the folder structures and naming convention might change requiring a revisit here
                            unpack_archive(
                                f"{output_folder}/Results/u2018_clc2012_v2020_20u1_raster100m.zip",
                                f"{output_folder}/Results/",
                            )
                            copy2(
                                f"{output_folder}/Results/u2018_clc2012_v2020_20u1_raster100m/DATA/U2018_CLC2012_V2020_20u1.tif",
                                tif_file,
                            )
                            break

                        else:
                            logger.info(f"Access error {file_response.status_code}")
            break
