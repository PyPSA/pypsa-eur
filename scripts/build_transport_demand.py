# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build land transport demand per clustered model region including efficiency
improvements due to drivetrain changes, time series for electric vehicle
availability and demand-side management constraints.
"""

import logging

import numpy as np
import pandas as pd
import xarray as xr
from _helpers import (
    configure_logging,
    generate_periodic_profiles,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

transport_cols = {"light":['Powered two-wheelers', 'Passenger cars',
                   'Light commercial vehicles'],
                  "heavy": ['Motor coaches, buses and trolley buses',
                                     'Heavy goods vehicles']}


def build_nodal_transport_data(fn, pop_layout, year):
    # get numbers of car and fuel efficiency per country
    transport_data = pd.read_csv(fn, index_col=[0, 1])
    transport_data = transport_data.xs(year, level="year")

    nodal_transport_data = transport_data.loc[pop_layout.ct].fillna(0.0)
    nodal_transport_data.index = pop_layout.index
    
    eff_cols = transport_data.columns.str.contains("efficiency")
    car_cols = transport_data.columns[~eff_cols]
    nodal_transport_data[car_cols] = (
       nodal_transport_data[car_cols].mul(pop_layout["fraction"], axis=0)
    )
    # fill missing fuel efficiency [kWh/100 km] with average data
    for col in nodal_transport_data.columns[eff_cols]:
        nodal_transport_data.loc[
            nodal_transport_data[col] == 0.0,
            col,
        ] = transport_data[col].mean()
    
    return nodal_transport_data


def get_shape(traffic_fn):
    traffic = pd.read_csv(traffic_fn, skiprows=2, usecols=["count"]).squeeze("columns")

    # create annual profile take account time zone + summer time
    transport_shape = generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=traffic.values,
    )
    transport_shape = transport_shape / transport_shape.sum()
    
    return transport_shape

def build_transport_demand(traffic_fn_Pkw, traffic_fn_Lkw,
                           airtemp_fn, nodes, nodal_transport_data):
    """
    Returns transport demand per bus in unit km driven [100 km].
    """
    transport_shape_light =  get_shape(traffic_fn_Pkw)
    
    transport_shape_heavy =  get_shape(traffic_fn_Lkw)

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
    
    # adjust for the heating/cooling demand from ICE totals
    ice_correction_light = ((transport_shape_light * (1 + dd_ICE)).sum()
                            / transport_shape_light.sum())
    ice_correction_heavy = ((transport_shape_light * (1 + dd_ICE)).sum()
                            / transport_shape_heavy.sum())
    
    # non-electrified rail
    non_elec_rail = (1 - (pop_weighted_energy_totals["electricity rail"]
                          / pop_weighted_energy_totals["total rail"]))
    
    # total demand of driven vehicle-km [mio km]
    light_duty = nodal_transport_data[[f'mio km-driven {car_type}'
                                       for car_type in transport_cols["light"]]].sum(axis=1)
    heavy_duty = nodal_transport_data[[f'mio km-driven {car_type}'
                                       for car_type in transport_cols["heavy"]]].sum(axis=1)
    heavy_duty = pd.concat([heavy_duty,
                             non_elec_rail * nodal_transport_data['mio km-driven Rail']],
                            axis=1).sum(axis=1)
    
    

    def get_demand(profile, total, nyears, ice_correction, name):
        """Returns from total demand [mio km], given profile and ICE correction
        demand time-series in unit [100 km]."""

        demand = ((profile.multiply(total) * 1e4 * nyears)
                  .divide(ice_correction))

        return pd.concat([demand], keys=[name], axis=1)

    demand_light = get_demand(transport_shape_light, light_duty, nyears,
                              ice_correction_light, name="light")
    demand_heavy = get_demand(transport_shape_heavy, heavy_duty, nyears,
                              ice_correction_heavy, name="heavy")

    return pd.concat([demand_light, demand_heavy], axis=1)


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
    # car count in typical week
    traffic = pd.read_csv(fn, skiprows=2, usecols=["count"]).squeeze("columns")
    # maximum share plugged-in availability for passenger electric vehicles
    avail_max = options["bev_avail_max"]
    # average share plugged-in availability for passenger electric vehicles
    avail_mean = options["bev_avail_mean"]

    # linear scaling, highest when traffic is lowest, decreases if traffic increases
    avail = avail_max - (avail_max - avail_mean) * (traffic - traffic.min()) / (
        traffic.mean() - traffic.min()
    )

    if not avail[avail < 0].empty:
        logger.warning(
            "The BEV availability weekly profile has negative values which can "
            "lead to infeasibility."
        )

    return generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=avail.values,
    )


def bev_dsm_profile(snapshots, nodes, options):
    dsm_week = np.zeros((24 * 7,))

    # assuming that at a certain time ("bev_dsm_restriction_time") EVs have to
    # be charged to a minimum value (defined in bev_dsm_restriction_value)
    dsm_week[(np.arange(0, 7, 1) * 24 + options["bev_dsm_restriction_time"])] = options[
        "bev_dsm_restriction_value"
    ]

    return generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=dsm_week,
    )


def build_registrations(nodal_transport_data):
    share_reg = {}
    for transport_type in ["light", "heavy"]:
        cols = transport_cols[transport_type]
        number = nodal_transport_data[[f"Number {car_type}" for car_type in cols]]
        new_reg = nodal_transport_data[[f"New registration {car_type}" for car_type in cols]]
        share_reg[transport_type] = new_reg.sum(axis=1).div(number.sum(axis=1))
   
    pd.concat(share_reg).to_csv(snakemake.output.car_registration)
    
# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_transport_demand",
            simpl="",
            clusters=37,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    nodes = pop_layout.index

    pop_weighted_energy_totals = pd.read_csv(
        snakemake.input.pop_weighted_energy_totals, index_col=0
    )

    options = snakemake.params.sector

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day, tz="UTC"
    )

    nyears = len(snapshots) / 8760

    energy_totals_year = snakemake.params.energy_totals_year
    nodal_transport_data = build_nodal_transport_data(
        snakemake.input.transport_data, pop_layout, energy_totals_year
    )

    transport_demand = build_transport_demand(
        snakemake.input.traffic_data_Pkw,
        snakemake.input.traffic_data_Lkw,
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
    
    build_registrations(nodal_transport_data)
