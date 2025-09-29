# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected helpers in scripts/add_electricity.py."""

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

from scripts.add_electricity import attach_load


def test_attach_load_uses_sanitized_busmap(tmp_path):
    """Loads are attached using trimmed busmap labels."""

    times = pd.date_range("2000-01-01", periods=2, freq="H")
    buses_raw = ["bus1", "bus2"]

    data = xr.DataArray(
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        coords={"time": times, "name": buses_raw},
        dims=["time", "name"],
        name="load",
    )

    load_path = tmp_path / "load.nc"
    data.to_netcdf(load_path)

    busmap_path = tmp_path / "busmap.csv"
    pd.DataFrame(
        {
            "Bus": ["bus1", "bus2"],
            "busmap": [" zone_1 ", "zone_2  "],
        }
    ).to_csv(busmap_path, index=False)

    n = pypsa.Network()
    n.set_snapshots(times)
    n.add("Bus", "zone_1")
    n.add("Bus", "zone_2")

    attach_load(n, load_path.as_posix(), busmap_path.as_posix(), scaling=1.0)

    assert sorted(n.loads.index) == ["zone_1", "zone_2"]
    assert list(n.loads_t.p_set.columns) == ["zone_1", "zone_2"]
