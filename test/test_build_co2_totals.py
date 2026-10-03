# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pytest

from scripts.build_co2_totals import build_co2_totals, ipcc_to_sector


@pytest.mark.parametrize(
    "code, sector",
    [
        ("1.A.1.a", "electricity"),
        ("1.A.1.bc", "industrial non-elec"),
        ("1.A.3.b", "road non-elec"),
        ("1.A.3.e", "industrial non-elec"),
        ("1.A.4", "buildings non-elec"),
        ("1.B.2", "industrial non-elec"),
        ("2.A.1", "industrial processes"),
        ("3.C.3", "agriculture"),
        ("4.C", "waste management"),
        ("5.B", "industrial non-elec"),
    ],
)
def test_ipcc_to_sector(code, sector):
    assert ipcc_to_sector(code) == sector


def test_ipcc_to_sector_unknown():
    with pytest.raises(ValueError):
        ipcc_to_sector("9.Z")


@pytest.fixture
def edgar_file(tmp_path):
    rows = [
        ("DEU", "1.A.1.a", 1000.0, 2000.0),
        ("DEU", "1.A.2", 500.0, None),
        ("DEU", "1.B.1", 500.0, 100.0),
        ("SCG", "1.A.1.a", 9316.0, 0.0),
        ("AIR", "1.A.3.a", 7000.0, 7000.0),
    ]
    df = pd.DataFrame(
        rows,
        columns=[
            "Country_code_A3",
            "ipcc_code_2006_for_standard_report",
            "Y_1990",
            "Y_1991",
        ],
    )
    fn = tmp_path / "edgar.xlsx"
    df.to_excel(fn, sheet_name="IPCC 2006", startrow=9, index=False)
    return fn


def test_build_co2_totals(edgar_file):
    co2 = build_co2_totals(edgar_file, ["DE", "RS", "ME", "XK"])

    assert co2.loc[("DE", 1990), "electricity"] == pytest.approx(1.0)
    assert co2.loc[("DE", 1990), "industrial non-elec"] == pytest.approx(1.0)
    assert co2.loc[("DE", 1991), "industrial non-elec"] == pytest.approx(0.1)
    split = co2.xs(1990, level="year").loc[["RS", "ME", "XK"], "electricity"]
    assert split.sum() == pytest.approx(9.316)
    assert split["RS"] > split["XK"] > split["ME"]


def test_build_co2_totals_missing_country(edgar_file):
    with pytest.raises(ValueError, match="FR"):
        build_co2_totals(edgar_file, ["DE", "FR"])
