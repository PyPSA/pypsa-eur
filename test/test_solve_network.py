# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import copy

import pandas as pd
import pypsa
import pytest

from scripts.solve_network import (
    MATH,
    MATH_WRITERS,
    MathWriter,
    carbon,
    emission_stores,
    flexible_egs,
    register_math,
)

CONSTRAINTS = {
    "lossy_bidirectional_links": "Link_bidirectional_sync",
    "pipe_retrofit": "Link_pipe_retrofit",
    "carbon": "GlobalConstraint-CO2Limit",
    "safe": "GlobalConstraint-SAFE",
    "ccl": "GlobalConstraint-CCL min DE solar",
    "equity": "equity_min",
    "solar_potential": "Bus_solar_potential",
    "operational_reserve": "operational_reserve_margin",
    "gas_boiler_retrofit": "Link_gas_boiler_retrofit",
    "flexible_egs": "Bus_geothermal_reservoir_capacity",
    "import_limit": "GlobalConstraint-import_limit",
}


def register(n: pypsa.Network, writer: MathWriter, config: dict) -> None:
    assert writer(n, config, "2030")
    n.optimize.overrides.register(MATH / f"{writer.__name__}.yaml")


@pytest.fixture
def math_config(config: dict) -> dict:
    config = copy.deepcopy(config)
    config["solving"]["constraints"].update(SAFE=True, CCL=True, EQ="EQ0.7c")
    config["electricity"]["SAFE_reservemargin"] = 0.1
    config["electricity"]["operational_reserve"]["activate"] = True
    config["sector"]["enhanced_geothermal"]["enable"] = True
    config["sector"]["imports"].update(enable=True, limit=10.0)
    return config


@pytest.fixture
def math_network() -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(range(2))
    n.set_investment_periods([2030])
    n.add("Carrier", "co2", co2_emissions=1.0)
    n.add("Bus", ["DE0", "DE1", "DE0 heat"], country="DE")
    n.add("Bus", "co2 atmosphere", carrier="co2")
    n.add("Store", "co2 atmosphere", bus="co2 atmosphere", e_nom=1e6)
    n.add(
        "GlobalConstraint",
        "CO2Limit",
        type="co2_atmosphere",
        carrier_attribute="co2_emissions",
        constant=1e6,
    )
    n.add("Load", "DE0", bus="DE0", p_set=[1.0, 2.0])
    n.add("Load", "DE0 urban heat", bus="DE0 heat", p_set=[2.0, 1.0])
    n.add(
        "Generator",
        ["DE0 solar", "DE0 solar-hsat", "DE0 OCGT", "DE0 import H2"],
        bus="DE0",
        carrier=["solar", "solar-hsat", "OCGT", "import H2"],
        p_nom=[0.0, 1.0, 0.0, 1.0],
        p_nom_max=10.0,
        p_nom_extendable=[True, True, True, False],
        marginal_cost=1.0,
    )
    n.add(
        "Link",
        ["DE0 urban gas boiler", "DE0 urban retrofitted H2 boiler"],
        bus0="DE0",
        bus1="DE0 heat",
        carrier=["urban gas boiler", "urban retrofitted H2 boiler"],
        p_nom=[5.0, 0.0],
        p_nom_extendable=[False, True],
    )
    n.add(
        "Link",
        ["gas pipeline DE0 <-> DE1", "H2 pipeline retrofitted DE0 <-> DE1"],
        bus0="DE0",
        bus1="DE1",
        carrier=["gas pipeline", "H2 pipeline retrofitted"],
        p_nom=[10.0, 0.0],
        p_nom_extendable=True,
    )
    n.add(
        "Link",
        ["DC DE0-DE1", "DC DE0-DE1-reversed"],
        bus0=["DE0", "DE1"],
        bus1=["DE1", "DE0"],
        carrier="DC",
        p_nom_extendable=True,
    )
    n.add(
        "Link",
        "DE0 geothermal well",
        bus0="DE0",
        bus1="DE0 heat",
        carrier="geothermal heat",
        p_nom_extendable=True,
    )
    n.add(
        "StorageUnit",
        "DE0 geothermal reservoir",
        bus="DE0 heat",
        carrier="geothermal heat",
        p_nom_extendable=True,
    )
    n.links["reversed"] = n.links.index.str.endswith("-reversed")
    return n


@pytest.mark.parametrize("writer", MATH_WRITERS, ids=lambda w: w.__name__)
def test_math_writer_adds_its_constraint(
    math_network: pypsa.Network, math_config: dict, writer: MathWriter
) -> None:
    register(math_network, writer, math_config)
    math_network.optimize.create_model(multi_investment_periods=True)
    constraint = math_network.model.constraints[CONSTRAINTS[writer.__name__]]
    assert (constraint.labels >= 0).any()


def test_register_math_replaces_stored_math(
    math_network: pypsa.Network, math_config: dict
) -> None:
    for _ in range(2):
        register_math(math_network, math_config, "2030")
    assert list(math_network.optimize.overrides) == list(CONSTRAINTS)
    math_network.optimize.create_model(multi_investment_periods=True)


def test_carbon_bounds_period_emissions_and_budget(math_config: dict) -> None:
    n = pypsa.Network()
    n.set_snapshots(range(2))
    n.set_investment_periods([2030, 2040])
    n.investment_period_weightings["years"] = 10.0
    n.add("Carrier", "co2", co2_emissions=1.0)
    n.add("Bus", ["electricity", "gas"])
    n.add("Bus", "co2 atmosphere", carrier="co2")
    n.add("Store", "co2 atmosphere", bus="co2 atmosphere", e_nom=1e3)
    n.add("Load", "load", bus="electricity", p_set=1.0)
    n.add("Generator", "wind", bus="electricity", p_nom=1.0, marginal_cost=10.0)
    n.add("Generator", "gas", bus="gas", p_nom=1.0, marginal_cost=1.0)
    n.add(
        "Link",
        "gas plant",
        bus0="gas",
        bus1="electricity",
        bus2="co2 atmosphere",
        efficiency2=1.0,
        p_nom=1.0,
    )
    n.add(
        "GlobalConstraint",
        ["co2 2040", "budget 2030"],
        type=["co2_atmosphere", "Co2Budget"],
        carrier_attribute="co2_emissions",
        investment_period=[2040, 2030],
        constant=[1.5, 5.0],
    )
    register(n, carbon, math_config)
    n.optimize(solver_name="highs", multi_investment_periods=True)

    level = n.stores_t.e["co2 atmosphere"].groupby(level="period").last()
    emissions = level.diff().fillna(level)
    assert emissions.to_dict() == pytest.approx({2030: 0.5, 2040: 1.5})


def test_carbon_rejects_several_carrier_attributes(
    math_network: pypsa.Network, math_config: dict
) -> None:
    math_network.add(
        "GlobalConstraint",
        "CH4Bound",
        type="co2_atmosphere",
        carrier_attribute="ch4_emissions",
        constant=100.0,
    )
    with pytest.raises(ValueError, match="several carrier attributes"):
        carbon(math_network, math_config, None)


@pytest.mark.parametrize(
    "attr,value", [("e_initial", 5.0), ("e_initial_per_period", True)]
)
def test_emission_stores_rejects_seeded_level(
    math_network: pypsa.Network, attr: str, value: float | bool
) -> None:
    math_network.stores.loc["co2 atmosphere", attr] = value
    with pytest.raises(ValueError, match="non-zero level"):
        emission_stores(math_network, pd.Index(["co2"]))


def test_emission_stores_selects_accumulating_stores(
    math_network: pypsa.Network,
) -> None:
    stores = emission_stores(math_network, pd.Index(["co2"]))
    assert list(stores.index) == ["co2 atmosphere"]
    math_network.stores.loc["co2 atmosphere", "e_cyclic"] = True
    assert emission_stores(math_network, pd.Index(["co2"])).empty


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


def test_flexible_egs_bounds_reservoir_by_wells_at_its_bus(
    egs_network: pypsa.Network, math_config: dict
) -> None:
    register(egs_network, flexible_egs, math_config)
    egs_network.optimize(solver_name="highs")

    p_nom_opt = egs_network.storage_units.p_nom_opt
    assert p_nom_opt.to_dict() == {"A reservoir": 130.0, "B reservoir": 5.0}
