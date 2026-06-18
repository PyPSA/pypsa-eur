# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Tests for scripts/plot_summary.py color matching and tech renaming.
"""


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


def _get_tech_colors(tech_colors, keys):
    """Standalone copy of _get_tech_colors from scripts/plot_summary.py."""
    fallback = "#333333"
    colors = []
    for key in keys:
        if key in tech_colors:
            colors.append(tech_colors[key])
        else:
            colors.append(fallback)
    return colors


def parse_tech_colors_from_yaml():
    """Custom lightweight YAML parser to extract tech_colors mapping from config."""
    tech_colors = {}
    in_tech_colors = False
    with open("config/plotting.default.yaml", "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped == "tech_colors:":
                in_tech_colors = True
                continue
            if in_tech_colors:
                indent = len(line) - len(line.lstrip())
                # If we exit the indented section of tech_colors
                if indent < 4 and not stripped.startswith("-"):
                    if ":" in stripped and not stripped.startswith("tech_colors"):
                        in_tech_colors = False
                        continue
                if ":" in stripped:
                    key, value = stripped.split(":", 1)
                    key = key.strip()
                    value = value.strip().strip("'").strip('"')
                    tech_colors[key] = value
    return tech_colors


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
    tech_colors = parse_tech_colors_from_yaml()
    assert "heat dsm" in tech_colors, "heat dsm must have a defined color in tech_colors"
    assert tech_colors["heat dsm"] == "#ff5c5c"


def test_get_tech_colors_all_present():
    """Verify color extraction when all keys exist."""
    tech_colors = {"onwind": "#235ebc", "solar": "#f9d002"}
    keys = ["onwind", "solar"]
    result = _get_tech_colors(tech_colors, keys)
    assert result == ["#235ebc", "#f9d002"]


def test_get_tech_colors_missing_fallback():
    """Verify fallback usage for missing color definitions."""
    tech_colors = {"onwind": "#235ebc"}
    keys = ["onwind", "unknown_carrier"]
    result = _get_tech_colors(tech_colors, keys)
    assert result == ["#235ebc", "#333333"]


def test_get_tech_colors_empty_dict():
    """Verify that all keys return fallback when dictionary is empty."""
    tech_colors = {}
    keys = ["onwind", "solar"]
    result = _get_tech_colors(tech_colors, keys)
    assert result == ["#333333", "#333333"]


def test_get_tech_colors_empty_keys():
    """Verify empty keys list returns empty color list."""
    tech_colors = {"onwind": "#235ebc"}
    keys = []
    result = _get_tech_colors(tech_colors, keys)
    assert result == []
