# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Copy used configuration files and important scripts for archiving.
"""

from pathlib import Path
from shutil import copy

import yaml

files = {
    "config/config.yaml": "config.yaml",
    "Snakefile": "Snakefile",
    "scripts/solve_network.py": "solve_network.py",
    "scripts/prepare_sector_network.py": "prepare_sector_network.py",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("copy_config")

    basepath = Path(f"results/{snakemake.params.RDIR}config/")

    for f, name in files.items():
        copy(f, basepath / name)

    with open(basepath / "config.snakemake.yaml", "w") as yaml_file:
        yaml.dump(
            snakemake.config,
            yaml_file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
