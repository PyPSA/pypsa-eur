# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected helpers in scripts/add_electricity.py."""

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

from scripts.add_electricity import attach_load, attach_storageunits


def test_attach_load(tmp_path):
    """Clustered demand is attached per bus with scaling applied."""

    times = pd.date_range("2000-01-01", periods=2, freq="h")
    buses = ["zone_1", "zone_2"]
    values = np.array([[1.0, 2.0], [3.0, 4.0]])

    data = xr.DataArray(
        values,
        coords={"time": times, "bus": buses},
        dims=["time", "bus"],
        name="electricity demand (MW)",
    )
    load_path = tmp_path / "electricity_demand.nc"
    data.to_netcdf(load_path)

    n = pypsa.Network()
    n.set_snapshots(times)
    n.add("Bus", buses)

    attach_load(n, load_path.as_posix(), scaling=2.0)

    assert sorted(n.loads.index) == buses
    assert sorted(n.loads_t.p_set.columns) == buses
    np.testing.assert_allclose(n.loads_t.p_set[buses].values, 2.0 * values)


def test_attach_storageunits_energy_basis():
    """A dispatched-basis `max_hours` sizes the store to sustain it, at no extra cost."""
    costs = pd.DataFrame(
        {
            "capital_cost": {"iron-air": 157959.0, "battery": 60630.0},
            "marginal_cost": {"iron-air": 0.0, "battery": 0.0},
            "lifetime": {"iron-air": 17.5, "battery": 17.5},
            "efficiency": {
                "iron-air battery charge": 0.74,
                "iron-air battery discharge": 0.63,
                "battery inverter": 0.96,
            },
        }
    )
    max_hours = {"iron-air": 100, "battery": 6}

    n = pypsa.Network()
    n.add("Bus", ["bus_1", "bus_2"])
    attach_storageunits(n, costs, n.buses.index, ["iron-air", "battery"], max_hours)

    su = n.storage_units.set_index("carrier")

    iron_air = su.loc["iron-air"]
    np.testing.assert_allclose(iron_air.max_hours, 158.73)
    # atol matches the 0.005 h the two-decimal rounding can cost.
    np.testing.assert_allclose(
        iron_air.max_hours * iron_air.efficiency_dispatch, 100, atol=0.005
    )
    np.testing.assert_allclose(iron_air.capital_cost, 157959.0)

    # Stored-basis carriers must be left alone, despite efficiency_dispatch < 1.
    np.testing.assert_allclose(su.loc["battery"].max_hours, 6)
