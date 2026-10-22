# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Tests for scripts/plot_summary.py color matching and tech renaming.
"""

import pytest

from scripts.plot_summary import check_tech_colors


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
