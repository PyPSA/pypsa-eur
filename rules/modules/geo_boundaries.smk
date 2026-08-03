# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

MODULE_NAME = "geo_boundaries"


module geo_boundaries:
    pathvars:
        shapes=resources("shapes/{scenario}.parquet"),
        logs=logs(MODULE_NAME),
        resources=resources(MODULE_NAME),
        results=resources(f"{MODULE_NAME}/results"),
    snakefile:
        github(
            "modelblocks-org/module_geo_boundaries",
            path="workflow/Snakefile",
            tag=config["modules"][MODULE_NAME]["version"],
        )
    config:
        config["modules"][MODULE_NAME]["config"]


use rule * from geo_boundaries exclude all as geo_boundaries_*
