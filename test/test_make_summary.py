# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pypsa
import pytest
from pypsa import NetworkCollection

from scripts.make_summary import (
    _loop_over_collection,
    assign_carriers,
    assign_locations,
    calculate_capacities,
    calculate_costs,
    calculate_cumulative_costs,
    calculate_energy,
    calculate_metrics,
    calculate_prices,
)


@pytest.fixture
def simple_network() -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=3, freq="h"))
    n.snapshot_weightings.loc[:, :] = 1.0

    n.add("Bus", "bus1", carrier="AC", location="DE")
    n.add("Bus", "bus2", carrier="AC", location="FR")

    n.add("Carrier", "solar")
    n.add("Carrier", "gas")

    n.add(
        "Generator",
        "solar1",
        bus="bus1",
        carrier="solar",
        p_nom=100,
        p_nom_opt=150,
        marginal_cost=0,
        capital_cost=1000,
    )
    n.add(
        "Generator",
        "gas1",
        bus="bus2",
        carrier="gas",
        p_nom=50,
        p_nom_opt=50,
        marginal_cost=50,
        capital_cost=500,
    )

    n.add(
        "Line", "line1", bus0="bus1", bus1="bus2", s_nom=100, s_nom_opt=120, length=100
    )

    n.generators_t.p = pd.DataFrame(
        {"solar1": [80, 100, 60], "gas1": [20, 0, 40]},
        index=n.snapshots,
    )
    n.buses_t.marginal_price = pd.DataFrame(
        {"bus1": [30, 20, 40], "bus2": [35, 25, 45]},
        index=n.snapshots,
    )

    return n


@pytest.fixture
def simple_network_collection(simple_network) -> NetworkCollection:
    n1 = simple_network.copy()
    n2 = simple_network.copy()
    n2.generators.p_nom_opt *= 1.5
    return NetworkCollection([n1, n2], index=pd.Index([2030, 2040], name="horizon"))


class TestLoopOverCollectionDecorator:
    def test_single_network_passthrough(self, simple_network):
        @_loop_over_collection
        def dummy_func(n: pypsa.Network) -> pd.Series:
            return pd.Series({"a": 1, "b": 2})

        result = dummy_func(simple_network)
        assert isinstance(result, pd.Series)
        assert list(result.index) == ["a", "b"]

    def test_collection_returns_dataframe(self, simple_network_collection):
        @_loop_over_collection
        def dummy_func(n: pypsa.Network) -> pd.Series:
            return pd.Series({"a": n.generators.p_nom_opt.sum(), "b": 2})

        result = dummy_func(simple_network_collection)
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == [2030, 2040]
        assert result.loc["a", 2040] == 1.5 * result.loc["a", 2030]

    def test_decorator_preserves_function_name(self):
        @_loop_over_collection
        def my_special_func(n: pypsa.Network) -> pd.Series:
            return pd.Series()

        assert my_special_func.__name__ == "my_special_func"


class TestAssignCarriers:
    def test_assigns_ac_carrier_to_lines(self, simple_network):
        del simple_network.lines["carrier"]
        assign_carriers(simple_network)
        assert "carrier" in simple_network.lines.columns
        assert simple_network.lines.carrier.iloc[0] == "AC"

    def test_preserves_existing_carrier(self, simple_network):
        simple_network.lines["carrier"] = "DC"
        assign_carriers(simple_network)
        assert simple_network.lines.carrier.iloc[0] == "DC"


class TestAssignLocations:
    def test_assigns_location_to_generators(self, simple_network):
        assign_locations(simple_network)
        assert "location" in simple_network.generators.columns
        assert simple_network.generators.loc["solar1", "location"] == "DE"
        assert simple_network.generators.loc["gas1", "location"] == "FR"

    def test_assigns_location_to_lines(self, simple_network):
        assign_locations(simple_network)
        assert "location" in simple_network.lines.columns


class TestCalculateCapacities:
    def test_single_network_returns_series(self, simple_network):
        result = calculate_capacities(simple_network)
        assert isinstance(result, pd.Series)
        assert result.sum() > 0

    def test_collection_returns_dataframe(self, simple_network_collection):
        result = calculate_capacities(simple_network_collection)
        assert isinstance(result, pd.DataFrame)
        assert 2030 in result.columns
        assert 2040 in result.columns


class TestCalculateCosts:
    def test_single_network_returns_series(self, simple_network):
        result = calculate_costs(simple_network)
        assert isinstance(result, pd.Series)
        assert "capital" in result.index.get_level_values("cost")
        assert "marginal" in result.index.get_level_values("cost")

    def test_collection_returns_dataframe(self, simple_network_collection):
        result = calculate_costs(simple_network_collection)
        assert isinstance(result, pd.DataFrame)
        assert 2030 in result.columns
        assert 2040 in result.columns


class TestCalculateEnergy:
    def test_single_network_returns_series(self, simple_network):
        result = calculate_energy(simple_network)
        assert isinstance(result, pd.Series)

    def test_collection_returns_dataframe(self, simple_network_collection):
        result = calculate_energy(simple_network_collection)
        assert isinstance(result, pd.DataFrame)


class TestCalculateMetrics:
    def test_single_network_returns_series(self, simple_network):
        result = calculate_metrics(simple_network)
        assert isinstance(result, pd.Series)
        assert "line_volume_AC" in result.index
        assert "total costs" in result.index

    def test_collection_returns_dataframe(self, simple_network_collection):
        result = calculate_metrics(simple_network_collection)
        assert isinstance(result, pd.DataFrame)


class TestCalculatePrices:
    def test_single_network_returns_series(self, simple_network):
        result = calculate_prices(simple_network)
        assert isinstance(result, pd.Series)
        assert "AC" in result.index

    def test_collection_returns_dataframe(self, simple_network_collection):
        result = calculate_prices(simple_network_collection)
        assert isinstance(result, pd.DataFrame)


class TestCalculateCumulativeCosts:
    def test_single_horizon_returns_empty(self):
        costs_df = pd.DataFrame({2030: [100, 200]})
        horizons = pd.Index([2030])
        result = calculate_cumulative_costs(costs_df, horizons)
        assert result.empty

    def test_multiple_horizons_returns_discounted_costs(self):
        costs_df = pd.DataFrame({2030: [100], 2040: [100], 2050: [100]})
        horizons = pd.Index([2030, 2040, 2050])
        result = calculate_cumulative_costs(costs_df, horizons)

        assert "social_discount_rate" == result.index.name
        assert 2030 in result.columns
        assert 2040 in result.columns
        assert 2050 in result.columns
        assert "cumulative_cost" in result.columns

    def test_zero_discount_rate_no_discounting(self):
        costs_df = pd.DataFrame({2030: [100], 2040: [100]})
        horizons = pd.Index([2030, 2040])
        result = calculate_cumulative_costs(costs_df, horizons)

        assert result.loc[0.0, 2030] == 100
        assert result.loc[0.0, 2040] == 100

    def test_positive_discount_rate_reduces_future_costs(self):
        costs_df = pd.DataFrame({2030: [100], 2040: [100]})
        horizons = pd.Index([2030, 2040])
        result = calculate_cumulative_costs(costs_df, horizons)

        assert result.loc[0.05, 2030] == 100
        assert result.loc[0.05, 2040] < 100
