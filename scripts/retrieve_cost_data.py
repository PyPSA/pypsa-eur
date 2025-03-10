# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve cost data from ``technology-data``.
"""

import logging
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_cost_data", year=2030)
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    version = snakemake.params.version
    if "/" in version:
        baseurl = f"https://raw.githubusercontent.com/{version}/outputs/"
    else:
        baseurl = f"https://raw.githubusercontent.com/PyPSA/technology-data/{version}/outputs/"
    filepath = Path(snakemake.output[0])
    url = baseurl + filepath.name

    print(url)

    to_fn = Path(rootpath) / filepath

    print(to_fn)

    logger.info(f"Downloading technology data from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, to_fn, disable=disable_progress)



############################## PyPSA-Eur
#
# Load the costs, load the modifications, and add them on-the-fly
#

    import pandas as pd

    costs_update = snakemake.params.costs_update


    if costs_update['enable']:

        logger.info(f'########## [PyPSA-Spain] <retrieve_cost_data.py> INFO: Including cost updates detailed in {costs_update["costs_update_file"]}')
        ########## "to_fn" contains the file, e.g. resources/costs_2030.csv
        df = pd.read_csv(to_fn)
        original_columns = df.columns.tolist()

        ########## Load the updates
        df_updates = pd.read_csv(costs_update['costs_update_file'])

        ########## Add updates
        # Combinamos los DataFrames usando merge en las columnas "technology" y "parameter"
        # El resto de columnas, si tienen valores distintos, las duplica con "_orig" y "_new". Es el caso de "value"
        merged_df = pd.merge(df, df_updates, on=['technology', 'parameter'], how='left', suffixes=('_orig', '_new')) 
        # Actualizamos los valores en "value" en df con los valores de df_updates donde "technology" y "parameter" son iguales
        merged_df['value'] = merged_df['value_new'].fillna(merged_df['value_orig'])
        # Eliminamos las columnas auxiliares
        merged_df.drop(['value_new', 'value_orig'], axis=1, inplace=True)

        ########## Arrange columns order 
        merged_df = merged_df[original_columns]

        ########## Export the updated file
        merged_df.to_csv(to_fn, index=False)



    logger.info(f"Technology data available at at {to_fn}")
