"""
Retrieve gas infrastructure data from https://zenodo.org/record/4767098/files/IGGIELGN.zip
"""

import logging
from helper import progress_retrieve

import zipfile
from pathlib import Path

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('retrieve_gas_network_data')
        rootpath = '..'
    else:
        rootpath = '.'

    url = "https://zenodo.org/record/4767098/files/IGGIELGN.zip"

    # Save locations
    zip_fn = Path(f"{rootpath}/IGGIELGN.zip")
    to_fn = Path(f"{rootpath}/data/gas_network/scigrid-gas")

    logger.info(f"Downloading databundle from '{url}'.")
    progress_retrieve(url, zip_fn)

    logger.info(f"Extracting databundle.")
    zipfile.ZipFile(zip_fn).extractall(to_fn)

    zip_fn.unlink()

    logger.info(f"Gas infrastructure data available in '{to_fn}'.")
