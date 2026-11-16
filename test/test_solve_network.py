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
    n.add("Carrier", "geothermal heat")
    n.add("Bus", ["EU", "A", "B"], carrier="geothermal heat")
    n.add(
        "Link",
        ["A well 0", "A well 1", "B well 0", "A well 0-2025"],
        bus0="EU",
        bus1=["A", "A", "B", "A"],
        carrier="geothermal heat",
        p_nom=[0.0, 0.0, 0.0, 100.0],
        p_nom_extendable=[True, True, True, False],
        p_nom_min=[10.0, 20.0, 5.0, 0.0],
        p_nom_max=[10.0, 20.0, 5.0, 100.0],
    )
    n.add(
        "StorageUnit",
        ["A reservoir", "B reservoir"],
        bus=["A", "B"],
        carrier="geothermal heat",
        p_nom_extendable=True,
        capital_cost=-1.0,
    )
    return n


def test_add_flexible_egs_constraint_bounds_reservoir_by_wells_at_its_bus(
    egs_network: pypsa.Network,
) -> None:
    egs_network.optimize(
        extra_functionality=lambda n, _: add_flexible_egs_constraint(n),
        solver_name="highs",
    )

    p_nom_opt = egs_network.storage_units.p_nom_opt
    assert p_nom_opt.to_dict() == {"A reservoir": 130.0, "B reservoir": 5.0}
