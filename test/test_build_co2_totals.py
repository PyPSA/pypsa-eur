# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pytest

from scripts.build_co2_totals import POPULATION, build_co2_totals


def eea_rows(country, rows):
    return [
        dict(
            Country_code=country,
            Pollutant_name="CO2",
            Sector_code=c,
            Year=y,
            emissions=v,
        )
        for c, y, v in rows
    ]


@pytest.fixture
def inventories(tmp_path):
    eea = eea_rows(
        "DE",
        [
            ("1", 1990, 1000.0),
            ("1.A.1.a", 1990, 400.0),
            ("1.A.4.c", 1990, 100.0),
            ("1.D.1.b", 1990, 50.0),
            ("3", 1990, 20.0),
            ("1", "1985-1987", 9999.0),
        ],
    )
    enc = (
        eea_rows("AL", [("1", 1990, 2000.0)])
        + eea_rows("RS", [("1", 1990, 8000.0), ("1.A.1.a", 1990, 8000.0)])
        + eea_rows("ME", [("1", 1990, 0.0)])
    )
    uk = pd.DataFrame(
        {
            "GHG": ["CO2", "CO2", "CO2", "CH4"],
            "CRT category": ["1A1ai", "1A2a", "1D1a", "1A1ai"],
            "Year": [1990] * 4,
            "Emissions (MtCO2e)": [2.0, 1.0, 0.5, 9.0],
        }
    )
    fns = tmp_path / "eea.csv", tmp_path / "enc.csv", tmp_path / "uk.xlsx"
    pd.DataFrame(eea).to_csv(fns[0], index=False)
    pd.DataFrame(enc).to_csv(fns[1], index=False)
    uk.to_excel(fns[2], sheet_name="UK_by_source", index=False)
    return fns


def test_build_co2_totals(inventories):
    co2 = build_co2_totals(["DE", "GB", "RS", "XK", "BA"], *inventories, "CO2")
    de, gb = co2.loc[("DE", 1990)], co2.loc[("GB", 1990)]

    assert co2.index.unique("year").tolist() == [1990]
    # bunkers are not part of total energy and must not reduce industry
    assert de["industrial non-elec"] == pytest.approx(0.5)
    assert de["international navigation"] == pytest.approx(0.05)
    assert de["agriculture"] == pytest.approx(0.12)
    assert gb["electricity"] == pytest.approx(2.0)
    assert gb["industrial non-elec"] == pytest.approx(1.0)
    assert gb["international aviation"] == pytest.approx(0.5)

    rs = co2.xs(1990, level="year").loc[["RS", "XK"], "electricity"]
    assert rs.sum() == pytest.approx(8.0)
    assert rs["RS"] / rs["XK"] == pytest.approx(POPULATION["RS"] / POPULATION["XK"])
    balkans = sum(POPULATION[c] for c in ["AL", "ME", "RS", "XK"])
    ba = co2.loc[("BA", 1990)]
    assert ba["electricity"] == pytest.approx(8.0 * POPULATION["BA"] / balkans)
    assert ba["industrial non-elec"] == pytest.approx(2.0 * POPULATION["BA"] / balkans)


def test_build_co2_totals_missing_country(inventories):
    with pytest.raises(ValueError, match="FR"):
        build_co2_totals(["DE", "FR"], *inventories, "CO2")
