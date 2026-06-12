# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pypsa
import pytest

from scripts._helpers import rename_network_component


def _build_network_with_generators() -> pypsa.Network:
    network = pypsa.Network()
    network.add("Bus", "bus")
    network.snapshots = pd.RangeIndex(2)
    network.add("Generator", "G1", bus="bus", p_nom=10)
    network.add("Generator", "G2", bus="bus", p_nom=20)
    network.generators_t.p = pd.DataFrame(
        [[0.5, 1.0], [0.7, 1.2]],
        index=network.snapshots,
        columns=["G1", "G2"],
    )
    return network


def test_rename_network_component_updates_static_and_dynamic_tables():
    network = _build_network_with_generators()
    rename_map = pd.Series(["G1_new", "G2_new"], index=["G1", "G2"])

    rename_network_component(network, "Generator", rename_map)

    static_index = network.static("Generator").index
    assert "G1_new" in static_index and "G1" not in static_index
    assert "G2_new" in static_index and "G2" not in static_index

    generator_p = network.generators_t.p
    assert list(generator_p.columns) == ["G1_new", "G2_new"]

    dynamic_tables = network.dynamic("Generator")
    assert generator_p is dynamic_tables["p"]


def test_rename_network_component_detects_name_conflicts():
    network = _build_network_with_generators()
    rename_map = pd.Series(["G2"], index=["G1"])

    with pytest.raises(ValueError, match="collide"):
        rename_network_component(network, "Generator", rename_map)


def test_rename_network_component_rejects_missing_keys():
    network = _build_network_with_generators()
    rename_map = pd.Series(["G3_new"], index=["missing"])

    with pytest.raises(KeyError, match="missing keys"):
        rename_network_component(network, "Generator", rename_map)
