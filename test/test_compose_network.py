# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected functions in scripts/compose_network.py."""

import pandas as pd
import pypsa

from scripts.compose_network import adjust_renewable_capacity_limits

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
