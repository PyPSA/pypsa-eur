# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve TABULA building typology data used in :mod:`build_retro_cost`.

The source is the TABULA calculator workbook published by EPISCOPE.
The ``Calc.Set.Building`` sheet is exported to the CSV format historically
shipped in the PyPSA-Eur data bundle.

When an archive mirror is added on data.pypsa.org, it should store the same
Excel workbook; both primary and archive sources will export the sheet to CSV.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

SHEET_NAME = "Calc.Set.Building"


def export_calcsetbuilding_csv(xlsx: Path, csv: Path) -> None:
    """Export the Calc.Set.Building sheet from the TABULA workbook to CSV.
    
    Parameters
    ----------
    xlsx : Path
        Path to the TABULA workbook.
    csv : Path
        Path to the CSV file to export.

    Returns
    -------
    None
    """
    logger.info("Exporting TABULA sheet '%s' from '%s'.", SHEET_NAME, xlsx)
    data = pd.read_excel(xlsx, sheet_name=SHEET_NAME, header=None)
    csv.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(csv, index=False, header=False)
    logger.info("Wrote '%s' (%s rows, %s columns).", csv, *data.shape)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tabula_calcsetbuilding")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    input_path = Path(snakemake.input.data)
    output_path = Path(snakemake.output.csv)

    export_calcsetbuilding_csv(input_path, output_path)
