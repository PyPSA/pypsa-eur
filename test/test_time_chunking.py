# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import pytest

from scripts._helpers import (
    build_time_chunked_snapshot_weightings,
    parse_time_chunking_resolution,
    update_config_from_wildcards,
)


class Wildcards(dict):
    __getattr__ = dict.__getitem__


def test_time_chunking_selects_equally_spaced_chunks_and_preserves_weighting():
    snapshots = pd.date_range("2013-01-01", periods=8760, freq="h")
    snapshot_weightings = pd.DataFrame(
        1.0, index=snapshots, columns=["objective", "generators", "stores"]
    )

    chunked = build_time_chunked_snapshot_weightings(
        snapshots, snapshot_weightings, "12c24"
    )

    assert len(chunked) == 12 * 24
    assert chunked.index.is_monotonic_increasing
    assert chunked.index[:24].equals(snapshots[:24])
    assert chunked.index[24] == snapshots[730]
    assert chunked.sum().eq(8760).all()
    assert chunked.objective.unique()[0] == pytest.approx(8760 / (12 * 24))


def test_time_chunking_wildcard_updates_temporal_resolution():
    config = {
        "clustering": {
            "temporal": {"resolution_elec": False, "resolution_sector": False}
        },
        "electricity": {"co2base": 1},
        "costs": {"emission_prices": {}},
        "autarky": {},
        "adjustments": {"electricity": False, "sector": False},
        "sector": {},
        "solving": {"constraints": {}},
        "lines": {},
        "links": {},
    }

    update_config_from_wildcards(config, Wildcards(opts="12c24", sector_opts=""))
    update_config_from_wildcards(config, Wildcards(opts="", sector_opts="12c24"))

    assert parse_time_chunking_resolution("12c24") == (12, 24)
    assert config["clustering"]["temporal"]["resolution_elec"] == "12c24"
    assert config["clustering"]["temporal"]["resolution_sector"] == "12c24"

