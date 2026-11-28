# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for the composites in scripts/composites."""

import pandas as pd
import pypsa
import pytest

from scripts.add_electricity import attach_stores
from scripts.prepare_sector_network import add_tes

COSTS = pd.DataFrame(
    {
        "capital_cost": {"battery storage": 10.0, "battery inverter": 10.0},
        "marginal_cost": 0.0,
        "lifetime": 20.0,
        "efficiency": {"battery storage": 1.0, "battery inverter": 0.81},
    }
)


@pytest.fixture
def n() -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(range(4))
    n.add("Bus", ["north", "north heat"])
    n.add("Generator", "sun", bus="north", p_nom=100, p_max_pu=[1, 1, 0, 0])
    n.add("Generator", "gas", bus="north", p_nom=100, marginal_cost=100)
    n.add("Link", "heat pump", bus0="north", bus1="north heat", p_nom=100)
    n.add("Load", ["power", "heat"], bus=["north", "north heat"], p_set=10)
    return n


def test_battery_charger_is_rated_by_its_discharger(n):
    attach_stores(n, COSTS, pd.Index(["north"]), ["battery"])
    n.optimize(solver_name="highs")
    links = n.c.links.static
    assert {"north battery charger", "north battery discharger"} <= set(links.index)
    assert links.carrier["north battery charger"] == "battery charger"
    p_nom = links.p_nom_opt
    assert p_nom["north battery discharger"] > 0
    assert p_nom["north battery charger"] == pytest.approx(
        0.9 * p_nom["north battery discharger"]
    )


def test_tes_ties_its_store_and_links(n):
    add_tes(
        n,
        pd.Index(["north"]),
        "water tanks",
        bus="north heat",
        discharger_efficiency=0.8,
        charger_ratio=0.8,
        store_capital_cost=1.0,
        lifetime=20.0,
        energy_to_power_ratio=6.0,
    )
    n.optimize(solver_name="highs")
    store = n.c.stores.static.loc["north water tanks"]
    p_nom = n.c.links.static.p_nom_opt
    assert store.carrier == "water tanks"
    assert p_nom["north water tanks discharger"] > 0
    assert p_nom["north water tanks charger"] == pytest.approx(
        0.8 * p_nom["north water tanks discharger"]
    )
    assert store.e_nom_opt == pytest.approx(6 * p_nom["north water tanks charger"])


def test_fixed_members_are_not_tied(n):
    attach_stores(n, COSTS, pd.Index(["north"]), ["battery"])
    members = n.composites["battery_inverter"].members.loc["north battery"]
    n.c.links.static.loc[members, ["p_nom", "p_nom_extendable"]] = [5.0, False]
    n.optimize(solver_name="highs")
    assert (n.c.links.static.p_nom_opt[members] == 5.0).all()
