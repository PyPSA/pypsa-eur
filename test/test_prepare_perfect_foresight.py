# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for perfect foresight utility functions in scripts/prepare_perfect_foresight.py."""

import numpy as np
import pandas as pd
import pypsa
import pytest

from scripts.prepare_perfect_foresight import (
    apply_investment_period_weightings,
    get_investment_weighting,
)


class TestGetInvestmentWeighting:
    """Tests for get_investment_weighting function."""

    def test_equal_periods_basic(self):
        """Test with three equal 10-year periods."""
        time_weighting = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
        result = get_investment_weighting(time_weighting, r=0.01)

        # Check that we get three weightings
        assert len(result) == 3
        # Check that earlier periods have higher weights (due to discounting)
        assert result[2030] > result[2040] > result[2050]
        # All values should be positive
        assert (result > 0).all()

    def test_equal_periods_zero_discount(self):
        """Test with zero discount rate - all periods should have equal weight."""
        time_weighting = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
        result = get_investment_weighting(time_weighting, r=0.0)

        # With zero discount, all periods have equal weight (equal to duration)
        np.testing.assert_allclose(result.values, [10.0, 10.0, 10.0], rtol=1e-10)

    def test_unequal_periods(self):
        """Test with unequal period lengths."""
        time_weighting = pd.Series([10, 15, 5], index=[2030, 2045, 2060])
        result = get_investment_weighting(time_weighting, r=0.02)

        # Check dimensions
        assert len(result) == 3
        # Second period has longer duration (15 years) so higher weight despite being later
        # Third period is both shorter and later, so lowest weight
        assert result[2045] > result[2030] > result[2060]

    def test_high_discount_rate(self):
        """Test with high discount rate (5%)."""
        time_weighting = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
        result_high = get_investment_weighting(time_weighting, r=0.05)
        result_low = get_investment_weighting(time_weighting, r=0.01)

        # Higher discount rate should make earlier periods more valuable
        # (ratio of 2030/2050 should be larger with higher discount)
        ratio_high = result_high[2030] / result_high[2050]
        ratio_low = result_low[2030] / result_low[2050]
        assert ratio_high > ratio_low

    def test_single_period(self):
        """Test with single period."""
        time_weighting = pd.Series([10], index=[2030])
        result = get_investment_weighting(time_weighting, r=0.02)

        assert len(result) == 1
        # Single period weight is discounted sum from year 0 to 10
        expected = sum(1 / (1.02) ** t for t in range(0, 10))
        np.testing.assert_allclose(result[2030], expected, rtol=1e-10)

    def test_two_periods(self):
        """Test with two periods."""
        time_weighting = pd.Series([10, 10], index=[2030, 2040])
        result = get_investment_weighting(time_weighting, r=0.02)

        assert len(result) == 2
        assert result[2030] > result[2040]

    def test_five_year_periods(self):
        """Test with shorter 5-year periods."""
        time_weighting = pd.Series([5, 5, 5, 5], index=[2030, 2035, 2040, 2045])
        result = get_investment_weighting(time_weighting, r=0.03)

        assert len(result) == 4
        # Should be monotonically decreasing
        assert all(result.iloc[i] > result.iloc[i + 1] for i in range(3))

    def test_typical_pypsa_eur_scenario(self):
        """Test with typical PyPSA-Eur planning horizons."""
        # Common scenario: 2030, 2040, 2050 with 10-year gaps
        time_weighting = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
        result = get_investment_weighting(time_weighting, r=0.01)

        # Verify result structure
        assert result.index.tolist() == [2030, 2040, 2050]
        # First period should have highest weight
        assert result[2030] == result.max()
        # Last period should have lowest weight
        assert result[2050] == result.min()

    def test_weightings_sum_property(self):
        """Test that weightings have reasonable sum (< total time horizon)."""
        time_weighting = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
        result = get_investment_weighting(time_weighting, r=0.02)

        # Sum of weighted years should be less than total years due to discounting
        total_years = time_weighting.sum()
        weighted_sum = result.sum()
        assert weighted_sum < total_years


class TestApplyInvestmentPeriodWeightings:
    """Tests for apply_investment_period_weightings function."""

    def test_basic_application(self):
        """Test basic application to multi-period network."""
        n = pypsa.Network()
        # Create multi-period network
        snapshots = pd.MultiIndex.from_product(
            [[2030, 2040, 2050], pd.date_range("2030-01-01", periods=3, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots)

        apply_investment_period_weightings(n, social_discountrate=0.02)

        # Check that weightings were added
        assert "years" in n.investment_period_weightings.columns
        assert "objective" in n.investment_period_weightings.columns

        # Check dimensions
        assert len(n.investment_period_weightings) == 3

        # Check years column (should be 10, 10, 10 with ffill for last)
        np.testing.assert_array_equal(
            n.investment_period_weightings["years"].values, [10.0, 10.0, 10.0]
        )

        # Check objective weightings decrease
        obj = n.investment_period_weightings["objective"]
        assert obj[2030] > obj[2040] > obj[2050]

    def test_raises_on_empty_investment_periods(self):
        """Test that function raises error for network without investment periods."""
        n = pypsa.Network()
        n.set_snapshots(pd.date_range("2030-01-01", periods=10, freq="h"))

        with pytest.raises(
            ValueError, match="Cannot apply investment period weightings"
        ):
            apply_investment_period_weightings(n, social_discountrate=0.02)

    def test_single_period_network(self):
        """Test with single investment period."""
        n = pypsa.Network()
        snapshots = pd.MultiIndex.from_product(
            [[2030], pd.date_range("2030-01-01", periods=5, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots)

        apply_investment_period_weightings(n, social_discountrate=0.02)

        # Should have weightings even for single period
        assert len(n.investment_period_weightings) == 1
        # Years should be filled (using ffill, so will be NaN converted to some value)
        assert not n.investment_period_weightings["years"].isna().all()

    def test_unequal_period_spacing(self):
        """Test with unevenly spaced investment periods."""
        n = pypsa.Network()
        snapshots = pd.MultiIndex.from_product(
            [[2030, 2045, 2060], pd.date_range("2030-01-01", periods=3, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots)

        apply_investment_period_weightings(n, social_discountrate=0.03)

        # Check years are calculated correctly
        years = n.investment_period_weightings["years"]
        assert years[2030] == 15.0  # 2045 - 2030
        assert years[2045] == 15.0  # 2060 - 2045
        assert years[2060] == 15.0  # ffill from previous

    def test_zero_discount_rate(self):
        """Test with zero discount rate."""
        n = pypsa.Network()
        snapshots = pd.MultiIndex.from_product(
            [[2030, 2040, 2050], pd.date_range("2030-01-01", periods=3, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots)

        apply_investment_period_weightings(n, social_discountrate=0.0)

        # With zero discount, objective weights should equal years
        obj = n.investment_period_weightings["objective"]
        years = n.investment_period_weightings["years"]
        np.testing.assert_allclose(obj.values, years.values, rtol=1e-10)

    def test_idempotency(self):
        """Test that applying weightings multiple times is idempotent."""
        n = pypsa.Network()
        snapshots = pd.MultiIndex.from_product(
            [[2030, 2040], pd.date_range("2030-01-01", periods=3, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots)

        apply_investment_period_weightings(n, social_discountrate=0.02)
        first_result = n.investment_period_weightings.copy()

        # Apply again
        apply_investment_period_weightings(n, social_discountrate=0.02)
        second_result = n.investment_period_weightings.copy()

        # Results should be identical
        pd.testing.assert_frame_equal(first_result, second_result)

    def test_different_discount_rates(self):
        """Test that different discount rates produce different weightings."""
        n1 = pypsa.Network()
        n2 = pypsa.Network()

        for n in [n1, n2]:
            snapshots = pd.MultiIndex.from_product(
                [[2030, 2040, 2050], pd.date_range("2030-01-01", periods=3, freq="h")],
                names=["period", "timestep"],
            )
            n.set_snapshots(snapshots)

        apply_investment_period_weightings(n1, social_discountrate=0.01)
        apply_investment_period_weightings(n2, social_discountrate=0.05)

        # Objective weightings should differ
        obj1 = n1.investment_period_weightings["objective"]
        obj2 = n2.investment_period_weightings["objective"]

        # They should not be equal
        assert not np.allclose(obj1.values, obj2.values)

        # Higher discount rate should favor earlier periods more
        ratio1 = obj1[2030] / obj1[2050]
        ratio2 = obj2[2030] / obj2[2050]
        assert ratio2 > ratio1

    def test_four_periods(self):
        """Test with four investment periods."""
        n = pypsa.Network()
        snapshots = pd.MultiIndex.from_product(
            [
                [2030, 2040, 2050, 2060],
                pd.date_range("2030-01-01", periods=3, freq="h"),
            ],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots)

        apply_investment_period_weightings(n, social_discountrate=0.02)

        assert len(n.investment_period_weightings) == 4
        obj = n.investment_period_weightings["objective"]
        # Should be strictly decreasing
        assert all(obj.iloc[i] > obj.iloc[i + 1] for i in range(3))

    def test_recalculation_on_growing_network(self):
        """Test that weightings change correctly when network grows (typical incremental use)."""
        # Simulate incremental composition: first 2030, then add 2040
        n = pypsa.Network()

        # First: single period [2030]
        snapshots1 = pd.MultiIndex.from_product(
            [[2030], pd.date_range("2030-01-01", periods=3, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots1)
        apply_investment_period_weightings(n, social_discountrate=0.02)
        weight_2030_alone = n.investment_period_weightings.loc[2030, "objective"]

        # Now add 2040 (simulate concatenation)
        snapshots2 = pd.MultiIndex.from_product(
            [[2030, 2040], pd.date_range("2030-01-01", periods=3, freq="h")],
            names=["period", "timestep"],
        )
        n.set_snapshots(snapshots2)
        apply_investment_period_weightings(n, social_discountrate=0.02)
        weight_2030_with_2040 = n.investment_period_weightings.loc[2030, "objective"]

        # Weight for 2030 should change when 2040 is added
        # (because T changes in the formula)
        assert weight_2030_alone != weight_2030_with_2040
