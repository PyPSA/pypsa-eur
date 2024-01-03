# -*- coding: utf-8 -*-
# Copyright 2023 Thomas Gilon (Climact)
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This rule downloads the existing capacities from `IRENASTAT <https://www.irena.org/Data/Downloads/IRENASTAT>`_ and extracts it in the ``data/existing_capacities`` sub-directory.

**Relevant Settings**

.. code:: yaml

    enable:
        retrieve_irena:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`enable_cf`

**Outputs**

- ``data/existing_capacities``: existing capacities for offwind, onwind and solar

"""

import logging

import pandas as pd
from _helpers import configure_logging

logger = logging.getLogger(__name__)

REGIONS = [
    "Albania",
    "Austria",
    "Belgium",
    "Bosnia and Herzegovina",
    "Bulgaria",
    "Croatia",
    "Czechia",
    "Denmark",
    "Estonia",
    "Finland",
    "France",
    "Germany",
    "Greece",
    "Hungary",
    "Ireland",
    "Italy",
    "Latvia",
    "Lithuania",
    "Luxembourg",
    "Montenegro",
    # "Netherlands",
    "Netherlands (Kingdom of the)",
    "North Macedonia",
    "Norway",
    "Poland",
    "Portugal",
    "Romania",
    "Serbia",
    "Slovakia",
    "Slovenia",
    "Spain",
    "Sweden",
    "Switzerland",
    # "United Kingdom",
    "United Kingdom of Great Britain and Northern Ireland (the)",
]

REGIONS_DICT = {
    "Bosnia and Herzegovina": "Bosnia Herzg",
    "Netherlands (Kingdom of the)": "Netherlands",
    "United Kingdom of Great Britain and Northern Ireland (the)": "UK",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_irena")
    configure_logging(snakemake)

    irena_raw = pd.read_csv(
        "https://pxweb.irena.org:443/sq/99e64b12-fe03-4a7b-92ea-a22cc3713b92",
        skiprows=2,
        index_col=[0, 1, 3],
        encoding="latin-1",
    )

    var = "Installed electricity capacity (MW)"
    irena = irena_raw[var].unstack(level=2).reset_index(level=1).replace(0, "")

    irena = irena[irena.index.isin(REGIONS)]
    irena.rename(index=REGIONS_DICT, inplace=True)

    df_offwind = irena[irena.Technology.str.contains("Offshore")].drop(
        columns=["Technology"]
    )
    df_onwind = irena[irena.Technology.str.contains("Onshore")].drop(
        columns=["Technology"]
    )
    df_pv = irena[irena.Technology.str.contains("Solar")].drop(columns=["Technology"])

    df_offwind.to_csv(snakemake.output["offwind"])
    df_onwind.to_csv(snakemake.output["onwind"])
    df_pv.to_csv(snakemake.output["solar"])
