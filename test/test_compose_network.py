# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected functions of the network composition pipeline."""

import pandas as pd
import pypsa

from scripts.add_brownfield import add_brownfield, adjust_renewable_capacity_limits
from scripts.prepare_perfect_foresight import concatenate_network_with_previous

# Default renewable carriers for testing
DEFAULT_RENEWABLE_CARRIERS = [
    "solar",
    "solar rooftop",
    "solar-hsat",
    "onwind",
    "offwind-ac",
    "offwind-dc",
    "offwind-float",
]


def test_adjust_renewable_capacity_limits_basic():
    """Test that existing capacities are correctly subtracted from p_nom_max."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    # Add buses
    n.add("Bus", "bus1")
    n.add("Bus", "bus2")

    # Add extendable generators for current horizon (2030)
    n.add(
        "Generator",
        "bus1 solar-2030",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=1000,
        p_nom_min=0,
    )
    n.add(
        "Generator",
        "bus2 onwind-2030",
        bus="bus2",
        carrier="onwind",
        p_nom_extendable=True,
        p_nom_max=800,
        p_nom_min=0,
    )

    # Add non-extendable generators representing existing capacity from previous horizon
    n.add(
        "Generator",
        "bus1 solar-2020",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=False,
        p_nom=200,
    )
    n.add(
        "Generator",
        "bus2 onwind-2025",
        bus="bus2",
        carrier="onwind",
        p_nom_extendable=False,
        p_nom=300,
    )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2030", ["solar", "onwind"])

    # Check that p_nom_max has been reduced by existing capacity
    assert n.generators.loc["bus1 solar-2030", "p_nom_max"] == 800  # 1000 - 200
    assert n.generators.loc["bus2 onwind-2030", "p_nom_max"] == 500  # 800 - 300


def test_adjust_renewable_capacity_limits_multiple_existing():
    """Test handling of multiple existing generators at the same bus."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2050-01-01", periods=2, freq="h"))

    # Add bus
    n.add("Bus", "bus1")

    # Add extendable generator for current horizon
    n.add(
        "Generator",
        "bus1 solar-2050",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=1000,
        p_nom_min=0,
    )

    # Add multiple non-extendable generators from different previous horizons
    n.add(
        "Generator",
        "bus1 solar-2020",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=False,
        p_nom=100,
    )
    n.add(
        "Generator",
        "bus1 solar-2030",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=False,
        p_nom=150,
    )
    n.add(
        "Generator",
        "bus1 solar-2040",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=False,
        p_nom=200,
    )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2050", ["solar"])

    # Check that all existing capacities are summed and subtracted
    expected = 1000 - (100 + 150 + 200)  # 550
    assert n.generators.loc["bus1 solar-2050", "p_nom_max"] == expected


def test_adjust_renewable_capacity_limits_exceeds_potential():
    """Test handling when existing capacity exceeds technical potential."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    # Add bus
    n.add("Bus", "bus1")

    # Add extendable generator with small p_nom_max
    n.add(
        "Generator",
        "bus1 solar-2030",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=500,
        p_nom_min=100,
    )

    # Add non-extendable generator with capacity exceeding p_nom_max
    n.add(
        "Generator",
        "bus1 solar-2020",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=False,
        p_nom=600,
    )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2030", ["solar"])

    # Check that p_nom_max is adjusted to p_nom_min when existing exceeds potential
    assert n.generators.loc["bus1 solar-2030", "p_nom_max"] == 100  # Set to p_nom_min


def test_adjust_renewable_capacity_limits_negative_clipping():
    """Test that negative p_nom_max values are clipped to zero."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    # Add bus
    n.add("Bus", "bus1")

    # Add extendable generator
    n.add(
        "Generator",
        "bus1 offwind-ac-2030",
        bus="bus1",
        carrier="offwind-ac",
        p_nom_extendable=True,
        p_nom_max=400,
        p_nom_min=0,
    )

    # Add non-extendable generator with larger capacity
    n.add(
        "Generator",
        "bus1 offwind-ac-2020",
        bus="bus1",
        carrier="offwind-ac",
        p_nom_extendable=False,
        p_nom=600,
    )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2030", ["offwind-ac"])

    # Check that p_nom_max is clipped to 0
    assert n.generators.loc["bus1 offwind-ac-2030", "p_nom_max"] == 0


def test_adjust_renewable_capacity_limits_all_carriers():
    """Test that all renewable carriers are handled correctly."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    carriers = [
        "solar",
        "solar rooftop",
        "solar-hsat",
        "onwind",
        "offwind-ac",
        "offwind-dc",
        "offwind-float",
    ]

    for i, carrier in enumerate(carriers):
        bus_name = f"bus{i}"
        n.add("Bus", bus_name)

        # Add extendable generator
        n.add(
            "Generator",
            f"{bus_name} {carrier}-2030",
            bus=bus_name,
            carrier=carrier,
            p_nom_extendable=True,
            p_nom_max=1000,
            p_nom_min=0,
        )

        # Add existing generator
        n.add(
            "Generator",
            f"{bus_name} {carrier}-2020",
            bus=bus_name,
            carrier=carrier,
            p_nom_extendable=False,
            p_nom=100,
        )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2030", DEFAULT_RENEWABLE_CARRIERS)

    # Check all carriers have reduced p_nom_max
    for i, carrier in enumerate(carriers):
        bus_name = f"bus{i}"
        gen_name = f"{bus_name} {carrier}-2030"
        assert n.generators.loc[gen_name, "p_nom_max"] == 900  # 1000 - 100


def test_adjust_renewable_capacity_limits_no_existing():
    """Test that function works correctly when there are no existing generators."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    # Add bus
    n.add("Bus", "bus1")

    # Add only extendable generators (no existing)
    n.add(
        "Generator",
        "bus1 solar-2030",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=1000,
        p_nom_min=0,
    )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2030", ["solar"])

    # Check that p_nom_max remains unchanged
    assert n.generators.loc["bus1 solar-2030", "p_nom_max"] == 1000


def test_adjust_renewable_capacity_limits_non_renewable_unaffected():
    """Test that non-renewable generators are not affected."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    # Add bus
    n.add("Bus", "bus1")

    # Add renewable generator
    n.add(
        "Generator",
        "bus1 solar-2030",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=1000,
        p_nom_min=0,
    )

    # Add non-renewable generator (e.g., gas)
    n.add(
        "Generator",
        "bus1 CCGT",
        bus="bus1",
        carrier="CCGT",
        p_nom_extendable=True,
        p_nom_max=500,
        p_nom_min=0,
    )

    # Apply the function
    adjust_renewable_capacity_limits(n, "2030", ["solar"])

    # Check that non-renewable generator is unaffected
    assert n.generators.loc["bus1 CCGT", "p_nom_max"] == 500


def test_adjust_renewable_capacity_limits_custom_carriers():
    """Test that only carriers in the config list are adjusted."""
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))

    # Add buses
    n.add("Bus", "bus1")
    n.add("Bus", "bus2")

    # Add extendable generators for different carriers
    n.add(
        "Generator",
        "bus1 solar-2030",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=True,
        p_nom_max=1000,
        p_nom_min=0,
    )
    n.add(
        "Generator",
        "bus2 onwind-2030",
        bus="bus2",
        carrier="onwind",
        p_nom_extendable=True,
        p_nom_max=800,
        p_nom_min=0,
    )

    # Add existing generators
    n.add(
        "Generator",
        "bus1 solar-2020",
        bus="bus1",
        carrier="solar",
        p_nom_extendable=False,
        p_nom=200,
    )
    n.add(
        "Generator",
        "bus2 onwind-2020",
        bus="bus2",
        carrier="onwind",
        p_nom_extendable=False,
        p_nom=300,
    )

    # Apply the function with only "solar" in the carriers list
    adjust_renewable_capacity_limits(n, "2030", ["solar"])

    # Check that only solar is adjusted, onwind is not
    assert n.generators.loc["bus1 solar-2030", "p_nom_max"] == 800  # 1000 - 200
    assert n.generators.loc["bus2 onwind-2030", "p_nom_max"] == 800  # unchanged


# ========== Tests for concatenate_network_with_previous ==========


def test_concatenate_basic_two_horizons():
    """Test basic concatenation of two single-period networks."""
    # Create first horizon network (2030)
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=24, freq="h"))
    n1.add("Bus", "bus1")
    n1.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030}
    n1.investment_periods = [2030]

    # Create second horizon network (2040)
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=24, freq="h"))
    n2.add("Bus", "bus1")
    n2.add("Bus", "bus2")
    n2.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n2.add("Generator", "gen2", bus="bus2", carrier="onwind", p_nom=200)
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate
    n_result = concatenate_network_with_previous(n1, n2, 2040)

    # Check investment periods
    assert hasattr(n_result, "investment_periods")
    assert list(n_result.investment_periods) == [2030, 2040]

    # Check snapshots are multi-indexed
    assert isinstance(n_result.snapshots, pd.MultiIndex)
    assert n_result.snapshots.names == ["period", "timestep"]
    assert len(n_result.snapshots) == 48  # 24 + 24

    # Check components
    assert len(n_result.generators) == 2
    assert "gen1" in n_result.generators.index
    assert "gen2" in n_result.generators.index
    assert len(n_result.buses) == 2


def test_concatenate_three_horizons():
    """Test concatenating three networks incrementally."""
    # First horizon (2030)
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=8, freq="h"))
    n1.add("Bus", "bus1")
    n1.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030}
    n1.investment_periods = [2030]

    # Second horizon (2040)
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=8, freq="h"))
    n2.add("Bus", "bus1")
    n2.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n2.add("Generator", "gen2", bus="bus1", carrier="onwind", p_nom=150)
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate first two
    n12 = concatenate_network_with_previous(n1, n2, 2040)

    # Third horizon (2050)
    n3 = pypsa.Network()
    n3.set_snapshots(pd.date_range("2050-01-01", periods=8, freq="h"))
    n3.add("Bus", "bus1")
    n3.add("Bus", "bus2")
    n3.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n3.add("Generator", "gen2", bus="bus1", carrier="onwind", p_nom=150)
    n3.add("Generator", "gen3", bus="bus2", carrier="offwind-ac", p_nom=300)
    n3.snapshot_weightings["objective"] = 1.0
    n3.investment_periods = [2050]

    # Concatenate all three
    n_result = concatenate_network_with_previous(n12, n3, 2050)

    # Check investment periods
    assert list(n_result.investment_periods) == [2030, 2040, 2050]

    # Check snapshots
    assert len(n_result.snapshots) == 24  # 8 + 8 + 8

    # Check components
    assert len(n_result.generators) == 3
    assert len(n_result.buses) == 2


def test_concatenate_with_time_series():
    """Test concatenation with time-varying attributes."""
    # First horizon with one generator
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=4, freq="h"))
    n1.add("Bus", "bus1")
    n1.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030}
    n1.investment_periods = [2030]

    # Second horizon with two generators
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=4, freq="h"))
    n2.add("Bus", "bus1")
    n2.add("Generator", "gen1", bus="bus1", carrier="solar", p_nom=100)
    n2.add("Generator", "gen2", bus="bus1", carrier="onwind", p_nom=150)
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate
    n_result = concatenate_network_with_previous(n1, n2, 2040)

    # Check that both generators exist
    assert "gen1" in n_result.generators.index
    assert "gen2" in n_result.generators.index
    assert len(n_result.generators) == 2

    # Check snapshots are properly concatenated
    assert len(n_result.snapshots) == 8  # 4 + 4 snapshots


def test_concatenate_with_stores():
    """Test concatenation includes stores from both networks."""
    # First horizon
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=4, freq="h"))
    n1.add("Bus", "bus1")
    n1.add("Store", "battery", bus="bus1", e_cyclic=True, e_nom=100)
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030}
    n1.investment_periods = [2030]

    # Second horizon with additional store
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=4, freq="h"))
    n2.add("Bus", "bus1")
    n2.add("Store", "battery", bus="bus1", e_cyclic=True, e_nom=100)
    n2.add(
        "Store",
        "co2_store",
        bus="bus1",
        carrier="co2 stored",
        e_cyclic=False,
        e_nom=1000,
    )
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate
    n_result = concatenate_network_with_previous(n1, n2, 2040)

    # Check that both stores exist
    assert "battery" in n_result.stores.index
    assert "co2_store" in n_result.stores.index
    assert len(n_result.stores) == 2


def test_concatenate_with_links():
    """Test concatenation with links including time-varying efficiency."""
    # First horizon
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=3, freq="h"))
    n1.add("Bus", "bus1")
    n1.add("Bus", "bus2")
    n1.add("Link", "link1", bus0="bus1", bus1="bus2", efficiency=0.9, p_nom=1000)
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030}
    n1.investment_periods = [2030]

    # Second horizon with new link
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=3, freq="h"))
    n2.add("Bus", "bus1")
    n2.add("Bus", "bus2")
    n2.add("Link", "link1", bus0="bus1", bus1="bus2", efficiency=0.9, p_nom=1000)
    n2.add("Link", "link2", bus0="bus2", bus1="bus1", efficiency=0.85, p_nom=800)
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate
    n_result = concatenate_network_with_previous(n1, n2, 2040)

    # Check links
    assert len(n_result.links) == 2
    assert "link1" in n_result.links.index
    assert "link2" in n_result.links.index


def test_concatenate_preserves_metadata():
    """Test that network metadata is preserved."""
    # First horizon
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=2, freq="h"))
    n1.add("Bus", "bus1")
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030, "custom_key": "value1"}
    n1.investment_periods = [2030]

    # Second horizon
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=2, freq="h"))
    n2.add("Bus", "bus1")
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate
    n_result = concatenate_network_with_previous(n1, n2, 2040)

    # Check that metadata from first network is preserved
    assert "custom_key" in n_result.meta


def test_concatenate_static_to_time_varying():
    """Test that static attributes are properly converted to time series when needed."""
    # First horizon with static marginal_cost
    n1 = pypsa.Network()
    n1.set_snapshots(pd.date_range("2030-01-01", periods=3, freq="h"))
    n1.add("Bus", "bus1")
    n1.add("Generator", "gen1", bus="bus1", carrier="gas", marginal_cost=50.0)
    n1.snapshot_weightings["objective"] = 1.0
    n1.meta = {"horizon": 2030}
    n1.investment_periods = [2030]

    # Second horizon with new gen2
    n2 = pypsa.Network()
    n2.set_snapshots(pd.date_range("2040-01-01", periods=3, freq="h"))
    n2.add("Bus", "bus1")
    n2.add("Generator", "gen1", bus="bus1", carrier="gas", marginal_cost=50.0)
    n2.add("Generator", "gen2", bus="bus1", carrier="solar", marginal_cost=0.0)
    n2.snapshot_weightings["objective"] = 1.0
    n2.investment_periods = [2040]

    # Concatenate
    n_result = concatenate_network_with_previous(n1, n2, 2040)

    # Check that gen2 was added with static marginal_cost
    assert "gen2" in n_result.generators.index
    assert n_result.generators.loc["gen2", "marginal_cost"] == 0.0


def test_add_brownfield_h2_retrofit_reduces_gas_pipeline():
    """Test that retrofitted H2 capacity is removed from the gas pipeline."""
    n = pypsa.Network()
    n.add("Bus", ["A", "B"])
    n.add("Link", "gas pipeline A-B", bus0="A", bus1="B", carrier="gas pipeline")
    n.links["p_nom"] = n.links["p_nom_max"] = 100.0
    n_p = n.copy()
    kwargs = dict(bus0="A", bus1="B", carrier="H2 pipeline retrofitted", lifetime=50)
    n_p.add("Link", "H2 pipeline retrofitted A-B-2030", build_year=2030, **kwargs)
    n_p.links["p_nom_opt"] = 30.0
    n.add("Link", "H2 pipeline retrofitted A-B-2040", build_year=2040, **kwargs)

    add_brownfield(n, n_p, 2040, h2_retrofit=True, h2_retrofit_capacity_per_ch4=0.6)

    assert n.links.at["gas pipeline A-B", "p_nom"] == 50.0
