# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

MODULE_NAME = "geo_boundaries"


module geo_boundaries:
    pathvars:
        shapes=f"resources/modules/{MODULE_NAME}/{{scenario}}.parquet",
        logs=f"logs/modules/{MODULE_NAME}",
        resources=f"data/modules/{MODULE_NAME}",
        results=f"resources/modules/{MODULE_NAME}/results",
    snakefile:
        github(
            "modelblocks-org/module_geo_boundaries",
            path="workflow/Snakefile",
            tag=config["modules"][MODULE_NAME]["version"],
        )
    config:
        module_config(MODULE_NAME)


use rule * from geo_boundaries exclude all as geo_boundaries_*
