# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Tests for individual NUTS country level configuration parsing in build_admin_shapes and busmap_for_admin_regions.
"""


def extract_country_level(admin_levels, countries):
    """
    Extracts individual country administrative levels.
    Mirrors the updated logic in build_admin_shapes and busmap_for_admin_regions.
    """
    return {
        k: v for k, v in admin_levels.get("countries", {}).items() if k in countries
    }


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


def test_regression_old_code_silently_drops_overrides():
    """
    Regression test for #2207: the OLD dict comprehension iterated over
    top-level admin_levels keys, silently producing an empty dict.
    """
    admin_levels = {"level": "bz", "countries": {"DE": 2}}
    countries = ["DE", "FR", "GB"]

    # OLD buggy code (what was there before):
    old_result = {
        k: v for k, v in admin_levels.items() if (k != "level") and (k in countries)
    }
    assert old_result == {}, "Sanity check: old code must produce empty dict"

    # NEW fixed code:
    new_result = extract_country_level(admin_levels, countries)
    assert new_result == {"DE": 2}, "Fix must extract per-country overrides"


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
