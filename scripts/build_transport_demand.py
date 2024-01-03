# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build land transport demand per clustered model region including efficiency
improvements due to drivetrain changes, time series for electric vehicle
availability and demand-side management constraints.
"""

import numpy as np
import pandas as pd
import xarray as xr
from _helpers import generate_periodic_profiles


def build_nodal_transport_data(fn, pop_layout):
    transport_data = pd.read_csv(fn, index_col=0)

    nodal_transport_data = transport_data.loc[pop_layout.ct].fillna(0.0)
    nodal_transport_data.index = pop_layout.index
    nodal_transport_data["number cars"] = (
        pop_layout["fraction"] * nodal_transport_data["number cars"]
    )
    nodal_transport_data.loc[
        nodal_transport_data["average fuel efficiency"] == 0.0,
        "average fuel efficiency",
    ] = transport_data["average fuel efficiency"].mean()

    return nodal_transport_data


def build_transport_demand(traffic_fn, airtemp_fn, nodes, nodal_transport_data):
    ## Get overall demand curve for all vehicles

    traffic = pd.read_csv(traffic_fn, skiprows=2, usecols=["count"]).squeeze("columns")

    transport_shape = generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=traffic.values,
    )
    transport_shape = transport_shape / transport_shape.sum()

    # electric motors are more efficient, so alter transport demand

    plug_to_wheels_eta = options["bev_plug_to_wheel_efficiency"]
    battery_to_wheels_eta = plug_to_wheels_eta * options["bev_charge_efficiency"]

    efficiency_gain = (
        nodal_transport_data["average fuel efficiency"] / battery_to_wheels_eta
    )

    # get heating demand for correction to demand time series
    temperature = xr.open_dataarray(airtemp_fn).to_pandas()

    # correction factors for vehicle heating
    dd_ICE = transport_degree_factor(
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["ICE_lower_degree_factor"],
        options["ICE_upper_degree_factor"],
    )

    dd_EV = transport_degree_factor(
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["EV_lower_degree_factor"],
        options["EV_upper_degree_factor"],
    )

    # divide out the heating/cooling demand from ICE totals
    # and multiply back in the heating/cooling demand for EVs
    ice_correction = (transport_shape * (1 + dd_ICE)).sum() / transport_shape.sum()

    energy_totals_transport = (
        pop_weighted_energy_totals["total road"]
        + pop_weighted_energy_totals["total rail"]
        - pop_weighted_energy_totals["electricity rail"]
    )

    return (
        (transport_shape.multiply(energy_totals_transport) * 1e6 * nyears)
        .divide(efficiency_gain * ice_correction)
        .multiply(1 + dd_EV)
    )


def transport_degree_factor(
    temperature,
    deadband_lower=15,
    deadband_upper=20,
    lower_degree_factor=0.5,
    upper_degree_factor=1.6,
):
    """
    Work out how much energy demand in vehicles increases due to heating and
    cooling.

    There is a deadband where there is no increase. Degree factors are %
    increase in demand compared to no heating/cooling fuel consumption.
    Returns per unit increase in demand for each place and time
    """

    dd = temperature.copy()

    dd[(temperature > deadband_lower) & (temperature < deadband_upper)] = 0.0

    dT_lower = deadband_lower - temperature[temperature < deadband_lower]
    dd[temperature < deadband_lower] = lower_degree_factor / 100 * dT_lower

    dT_upper = temperature[temperature > deadband_upper] - deadband_upper
    dd[temperature > deadband_upper] = upper_degree_factor / 100 * dT_upper

    return dd


def bev_availability_profile(fn, snapshots, nodes, options):
    """
    Derive plugged-in availability for passenger electric vehicles.
    """
    traffic = pd.read_csv(fn, skiprows=2, usecols=["count"]).squeeze("columns")

    avail_max = options["bev_avail_max"]
    avail_mean = options["bev_avail_mean"]

    avail = avail_max - (avail_max - avail_mean) * (traffic - traffic.min()) / (
        traffic.mean() - traffic.min()
    )

    return generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=avail.values,
    )


def bev_dsm_profile(snapshots, nodes, options):
    dsm_week = np.zeros((24 * 7,))

    dsm_week[(np.arange(0, 7, 1) * 24 + options["bev_dsm_restriction_time"])] = options[
        "bev_dsm_restriction_value"
    ]

    return generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=dsm_week,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_transport_demand",
            simpl="",
            clusters=48,
        )

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    nodes = pop_layout.index

    pop_weighted_energy_totals = pd.read_csv(
        snakemake.input.pop_weighted_energy_totals, index_col=0
    )

    options = snakemake.params.sector

    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots, tz="UTC")

    nyears = len(snapshots) / 8760

    nodal_transport_data = build_nodal_transport_data(
        snakemake.input.transport_data, pop_layout
    )

    transport_demand = build_transport_demand(
        snakemake.input.traffic_data_KFZ,
        snakemake.input.temp_air_total,
        nodes,
        nodal_transport_data,
    )

    avail_profile = bev_availability_profile(
        snakemake.input.traffic_data_Pkw, snapshots, nodes, options
    )

    dsm_profile = bev_dsm_profile(snapshots, nodes, options)

    nodal_transport_data.to_csv(snakemake.output.transport_data)
    transport_demand.to_csv(snakemake.output.transport_demand)
    avail_profile.to_csv(snakemake.output.avail_profile)
    dsm_profile.to_csv(snakemake.output.dsm_profile)
