# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Copy used configuration files and important scripts for archiving.
"""


import yaml

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("copy_config")

    with open(snakemake.output[0], "w") as yaml_file:
        yaml.dump(
            snakemake.config,
            yaml_file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
