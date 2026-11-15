# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pypsa
import pytest

from scripts.solve_network import (
    add_co2_atmosphere_constraint,
    add_flexible_egs_constraint,
    emission_stores,
)


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


@pytest.fixture(scope="function")
def egs_network() -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(pd.date_range("2013-01-01", periods=3, freq="h"))
    n.add("Carrier", "geothermal heat")
    n.add(
        "Bus",
        [
            "EU enhanced geothermal systems",
            "A geothermal heat surface",
            "B geothermal heat surface",
        ],
        carrier="geothermal heat",
    )
    n.add(
        "Link",
        [
            "A enhanced geothermal 0",
            "A enhanced geothermal 1",
            "B enhanced geothermal 0",
        ],
        bus0="EU enhanced geothermal systems",
        bus1=[
            "A geothermal heat surface",
            "A geothermal heat surface",
            "B geothermal heat surface",
        ],
        carrier="geothermal heat",
        p_nom_extendable=True,
        capital_cost=1.0,
    )
    n.add(
        "StorageUnit",
        ["A geothermal reservoir", "B geothermal reservoir"],
        bus=["A geothermal heat surface", "B geothermal heat surface"],
        carrier="geothermal heat",
        p_nom_extendable=True,
        capital_cost=1.0,
    )
    return n


def test_add_flexible_egs_constraint_bounds_reservoir_by_wells_at_its_bus(
    egs_network: pypsa.Network,
) -> None:
    egs_network.optimize.create_model()
    add_flexible_egs_constraint(egs_network)

    m = egs_network.model
    constraint = m.constraints["upper_bound_charging_capacity_of_geothermal_reservoir"]
    labels = pd.concat(
        [
            m["Link-p_nom"].labels.to_series(),
            m["StorageUnit-p_nom"].labels.to_series(),
        ]
    )

    def terms(bus: str) -> dict[str, float]:
        row = constraint.lhs.sel(bus=bus)
        names = pd.Series(labels.index, index=labels.values)
        return {
            names[v]: c for v, c in zip(row.vars.values, row.coeffs.values) if v != -1
        }

    assert list(constraint.indexes["bus"]) == [
        "A geothermal heat surface",
        "B geothermal heat surface",
    ]
    assert terms("A geothermal heat surface") == {
        "A geothermal reservoir": 1.0,
        "A enhanced geothermal 0": -1.0,
        "A enhanced geothermal 1": -1.0,
    }
    assert terms("B geothermal heat surface") == {
        "B geothermal reservoir": 1.0,
        "B enhanced geothermal 0": -1.0,
    }
    assert (constraint.sign == "<=").all()
    assert (constraint.rhs == 0.0).all()


def test_add_flexible_egs_constraint_skips_without_reservoirs(
    egs_network: pypsa.Network,
) -> None:
    egs_network.remove("StorageUnit", egs_network.storage_units.index)
    egs_network.optimize.create_model()
    add_flexible_egs_constraint(egs_network)

    assert (
        "upper_bound_charging_capacity_of_geothermal_reservoir"
        not in egs_network.model.constraints
    )
