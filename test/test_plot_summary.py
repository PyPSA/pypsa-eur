# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Tests for scripts/plot_summary.py color matching and tech renaming.
"""

import pytest
import yaml


def rename_techs(label: str) -> str:
    """Standalone copy of rename_techs from scripts/_helpers.py for dependency-free testing."""
    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral ",
    ]

    rename_if_contains = [
        "CHP",
        "gas boiler",
        "biogas",
        "solar thermal",
        "air heat pump",
        "ground heat pump",
        "resistive heater",
        "Fischer-Tropsch",
    ]

    rename_if_contains_dict = {
        "water tanks": "hot water storage",
        "retrofitting": "building retrofitting",
        "battery": "battery storage",
        "H2 for industry": "H2 for industry",
        "land transport fuel cell": "land transport fuel cell",
        "land transport oil": "land transport oil",
        "oil shipping": "shipping oil",
    }

    rename = {
        "solar": "solar PV",
        "Sabatier": "methanation",
        "offwind": "offshore wind",
        "offwind-ac": "offshore wind (AC)",
        "offwind-dc": "offshore wind (DC)",
        "offwind-float": "offshore wind (Float)",
        "onwind": "onshore wind",
        "ror": "hydroelectricity",
        "hydro": "hydroelectricity",
        "PHS": "hydroelectricity",
        "NH3": "ammonia",
        "co2 Store": "DAC",
        "co2 stored": "CO2 sequestration",
        "AC": "transmission lines",
        "DC": "transmission lines",
        "B2B": "transmission lines",
    }

    for ptr in prefix_to_remove:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old, new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old, new in rename.items():
        if old == label:
            label = new
    return label


def check_tech_colors(tech_colors, keys):
    """Standalone copy of check_tech_colors from scripts/plot_summary.py."""
    missing = [k for k in keys if k not in tech_colors]
    if missing:
        raise KeyError(
            f"The following technology carrier(s) do not have a defined color in the plotting configuration: {missing}"
        )


def test_regression_rename_techs_heat_dsm():
    """
    Regression test for #2108: Verify rename_techs maps different dsm
    options to 'heat dsm' and that 'heat dsm' has a color defined in the config.
    """
    # 1. Check rename mappings
    assert rename_techs("rural heat dsm") == "heat dsm"
    assert rename_techs("urban central heat dsm") == "heat dsm"
    assert rename_techs("residential urban decentral heat dsm") == "heat dsm"

    # 2. Check config has the key 'heat dsm'
    with open("config/plotting.default.yaml", "r") as f:
        config = yaml.safe_load(f)
    tech_colors = config["plotting"]["tech_colors"]

    assert "heat dsm" in tech_colors, "heat dsm must have a defined color in tech_colors"
    assert tech_colors["heat dsm"] == "#ff5c5c"


def test_check_tech_colors_all_present():
    """Verify that no KeyError is raised when all keys exist in tech_colors."""
    tech_colors = {"onwind": "#235ebc", "solar": "#f9d002"}
    keys = ["onwind", "solar"]
    check_tech_colors(tech_colors, keys)


def test_check_tech_colors_missing_raises():
    """Verify that check_tech_colors raises a KeyError detailing the missing keys."""
    tech_colors = {"onwind": "#235ebc"}
    keys = ["onwind", "unknown_carrier"]
    with pytest.raises(KeyError) as excinfo:
        check_tech_colors(tech_colors, keys)
    assert "unknown_carrier" in str(excinfo.value)
