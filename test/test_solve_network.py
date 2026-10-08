# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pypsa
import pytest

from scripts.solve_network import add_co2_atmosphere_constraint, emission_stores


@pytest.fixture(scope="function")
def co2_network() -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=3, freq="h"))
    n.add("Carrier", "co2", co2_emissions=1.0)
    n.add("Bus", "co2 atmosphere", carrier="co2")
    n.add("Store", "co2 atmosphere", bus="co2 atmosphere", e_nom=1e6, e_min_pu=-1)
    n.add("Bus", "electricity")
    n.add("Load", "load", bus="electricity", p_set=1.0)
    n.add("Generator", "gas", bus="electricity", p_nom=10.0, marginal_cost=1.0)
    return n


@pytest.mark.parametrize("sense", ["<=", ">=", "=="])
def test_add_co2_atmosphere_constraint_honours_sense(
    co2_network: pypsa.Network, sense: str
) -> None:
    co2_network.add(
        "GlobalConstraint",
        "CO2Bound",
        type="co2_atmosphere",
        carrier_attribute="co2_emissions",
        sense=sense,
        constant=100.0,
    )
    co2_network.optimize.create_model()
    add_co2_atmosphere_constraint(co2_network, co2_network.snapshots)

    constraint = co2_network.model.constraints["GlobalConstraint-CO2Bound"]
    expected = {"<=": "<=", ">=": ">=", "==": "="}[sense]
    assert (constraint.sign == expected).all()
    assert (constraint.rhs == 100.0).all()


def test_add_co2_atmosphere_constraint_rejects_unknown_sense(
    co2_network: pypsa.Network,
) -> None:
    co2_network.add(
        "GlobalConstraint",
        "CO2Bound",
        type="co2_atmosphere",
        carrier_attribute="co2_emissions",
        sense="<",
        constant=100.0,
    )
    co2_network.optimize.create_model()
    with pytest.raises(ValueError, match="Unsupported sense"):
        add_co2_atmosphere_constraint(co2_network, co2_network.snapshots)


@pytest.mark.parametrize(
    "attr,value", [("e_initial", 5.0), ("e_initial_per_period", True)]
)
def test_emission_stores_rejects_seeded_level(
    co2_network: pypsa.Network, attr: str, value: float | bool
) -> None:
    co2_network.stores.loc["co2 atmosphere", attr] = value
    with pytest.raises(ValueError, match="non-zero level"):
        emission_stores(co2_network, pd.Index(["co2"]))


def test_emission_stores_selects_accumulating_stores(
    co2_network: pypsa.Network,
) -> None:
    assert list(emission_stores(co2_network, pd.Index(["co2"])).index) == [
        "co2 atmosphere"
    ]
    co2_network.stores.loc["co2 atmosphere", "e_cyclic"] = True
    assert emission_stores(co2_network, pd.Index(["co2"])).empty
