# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Tests for individual NUTS country level configuration parsing in build_admin_shapes and busmap_for_admin_regions.
"""


from scripts._helpers import extract_country_level


def test_extract_country_level_with_overrides():
    """Verify that overrides are correctly extracted from nested countries dictionary."""
    admin_levels = {"level": "bz", "countries": {"DE": 2, "FR": 3, "IT": 1}}
    countries = ["DE", "FR", "GB"]
    result = extract_country_level(admin_levels, countries)
    assert result == {"DE": 2, "FR": 3}


def test_extract_country_level_without_countries_key():
    """Verify default behavior when 'countries' key is missing or empty."""
    admin_levels = {"level": "bz"}
    countries = ["DE", "FR"]
    result = extract_country_level(admin_levels, countries)
    assert result == {}


def test_extract_country_level_empty_countries_list():
    """Verify behavior when list of countries to keep is empty."""
    admin_levels = {"level": "bz", "countries": {"DE": 2}}
    countries = []
    result = extract_country_level(admin_levels, countries)
    assert result == {}


def test_level_key_does_not_leak():
    """Ensure top-level 'level' key never appears in country overrides."""
    admin_levels = {"level": 2, "countries": {"DE": 3}}
    countries = ["DE", "level"]  # adversarial: "level" in countries list
    result = extract_country_level(admin_levels, countries)
    assert "level" not in result  # must not pick up the top-level key
    assert result == {"DE": 3}


def test_countries_dict_explicitly_empty():
    """Empty countries dict in config returns empty result."""
    admin_levels = {"level": "bz", "countries": {}}
    countries = ["DE"]
    result = extract_country_level(admin_levels, countries)
    assert result == {}
