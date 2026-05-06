# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create profiles for road transport demand using measured data from vehicle monitoring by the German Federal Highway Research Institute (BASt).

This rule downloads the data files, extracts them, and then aggregates the data to weekly profiles for two vehicle types:
- "kfz": All motor vehicles (="Kraftfahrzeuge", i.e. cars, trucks, buses, motorcycles)
- "pkw": Passenger cars only (="Personenkraftwagen")

Outputs
-------

- ``data/mobility_profiles/build/<version>/kfz.csv``: Weekly profile for all motor vehicles (cars, trucks, buses, motorcycles).
- ``data/mobility_profiles/build/<version>/pkw.csv``: Weekly profile for passenger cars only.

**kfz.csv**

    ===================  ==========  ==========  =========================================================
    Field                Dimensions  Unit        Description
    ===================  ==========  ==========  =========================================================
    day                  day        day of week  Day of the week (0=Monday, 6=Sunday)
    -------------------  ----------  -----------  ---------------------------------------------------------
    hour                 hour       hour of day  Hour of the day (0-23)
    -------------------  ----------  -----------  ---------------------------------------------------------
    count                day, hour   --          Aggregated vehicle counts for all motor vehicles
                                                 (across all aggregated years and street types)
    -------------------  ----------  ----------  ---------------------------------------------------------
    n_counts             day, hour   --          Number of data points that were aggregated.
    ===================  ==========  ==========  =========================================================

**pkw.csv**
    ===================  ==========  ==========  =========================================================
    Field                Dimensions  Unit        Description
    ===================  ==========  ==========  =========================================================
    day                  day         day of week Day of the week (0=Monday, 6=Sunday)
    -------------------  ----------  ----------- ---------------------------------------------------------
    hour                 hour        hour of day Hour of the day (0-23)
    -------------------  ----------  ----------- ---------------------------------------------------------
    count                day, hour   --          Aggregated vehicle counts for passenger cars only
                                                 (across all aggregated years and street types)
    -------------------  ----------  ----------  ---------------------------------------------------------
    n_counts             day, hour   --          Number of data points that were aggregated.
    ===================  ==========  ==========  =========================================================

"""

import logging
from pathlib import Path
from shutil import unpack_archive

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_mobility_profiles")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Zip files have been downloaded by snakemake as storage objects
    # unzip them and move the extracted files into the `raw_files` folder
    zip_files = [Path(fp) for fp in snakemake.input["zip_files"]]
    output_folder = Path(snakemake.output["raw_files"])

    output_files = []
    for zip_file in zip_files:
        unpack_archive(zip_file, output_folder)
        output_file = output_folder / zip_file.name.replace(".zip", ".txt")
        output_files.append(output_file)
        logger.info(f"Unpacked {zip_file} to {output_file}")

    # We use this to add some information on the files used to create the output file later
    file_names = [f.name for f in zip_files]

    # Load each file as a dataframe and then concatenate them
    # (faster this way)
    output_files = sorted(output_files)
    dfs = []
    for output_file in output_files:
        logger.info(f"Processing {str(output_file.name)}")

        df = pd.read_csv(
            output_file,
            skiprows=0,
            delimiter=";",
            # Only load the columns we need:
            # 'Datum' is the date (YYMMDD)
            # 'Wotag' is the day of the week (1=Mon, 7=Sun)
            # 'Stunde' is the hour of the day (1-24)
            # 'KFZ|Pkw_R1' and 'KFZ|Pkw_R2' are the counts in each direction on the road
            #
            usecols=[
                "Datum",
                "Wotag",
                "Stunde",
                "KFZ_R1",
                "KFZ_R2",
                "Pkw_R1",
                "Pkw_R2",
            ],
        )

        dfs.append(df)
    vehicle_counts = pd.concat(dfs, ignore_index=True)

    # Turn 'Datum' into a datetime (YYYY-MM-DD) - we don't need it for the aggregation, but helpful for debugging/checks
    vehicle_counts["date"] = pd.to_datetime("20" + vehicle_counts["Datum"].astype(str))

    # Rename columns to English
    vehicle_counts = vehicle_counts.rename(
        columns={
            "Wotag": "day",
            "Stunde": "hour",
        }
    )

    # Aggregate data for both directions on the road
    vehicle_counts["kfz"] = vehicle_counts["KFZ_R1"] + vehicle_counts["KFZ_R2"]
    vehicle_counts["pkw"] = vehicle_counts["Pkw_R1"] + vehicle_counts["Pkw_R2"]

    # vehicles types to aggregate for and output files to write to
    vehicle_types = {
        "kfz": snakemake.output["kfz"],
        "pkw": snakemake.output["pkw"],
    }
    for vehicle_type, output_file in vehicle_types.items():
        logger.info(f"Aggregating and writing {vehicle_type} data to {output_file}")

        aggregated_data = vehicle_counts.groupby(["day", "hour"], as_index=False).agg(
            count=(vehicle_type, "sum"), n_counts=(vehicle_type, "count")
        )

        # Adjust day and hour to start from 0 (to match the expected data conventions)
        aggregated_data["day"] -= 1
        aggregated_data["hour"] -= 1

        aggregated_data.to_csv(output_file, index=False)

        # Add additional information to the beginning of the file
        # this information and format is similar to the original data format used by PyPSA-Eur previously
        with open(output_file, "r+") as file:
            content = file.read()
            file.seek(0, 0)
            file.write(
                f"# File generated for type: {vehicle_type} using data for: {', '.join(file_names)}\n"
            )
            file.write(f"# Time of generation: {pd.Timestamp.now()}\n")
            file.write(content)
