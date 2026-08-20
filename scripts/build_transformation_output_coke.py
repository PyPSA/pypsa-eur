# SPDX-FileCopyrightText: : 2025 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Extract transformation output for coke ovens from Eurostat energy balance data.

Outputs
-------
- ``resources/<run_name>/transformation_output_coke.csv``: Transformation output
  of coke ovens by country and energy carrier.
"""

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_transformation_output_coke")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    eurostat = pd.read_csv(snakemake.input.eurostat)
    eurostat.query("nrg_bal == 'TO_CO'").set_index(["country", "year", "siec"])[
        "value"
    ].unstack("siec").to_csv(snakemake.output.transformation_output_coke)
