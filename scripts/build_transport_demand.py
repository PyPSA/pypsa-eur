# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
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
import pypsa

from scripts._helpers import (
    configure_logging,
    generate_periodic_profiles,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def build_nodal_transport_data(fn, pop_layout, year):
    # get numbers of car and fuel efficiency per country
    transport_data = pd.read_csv(fn, index_col=[0, 1])
    transport_data = transport_data.xs(year, level="year")

    # break number of cars down to nodal level based on population density
    nodal_transport_data = transport_data.loc[pop_layout.ct].fillna(0.0)
    nodal_transport_data.index = pop_layout.index
    # add nodal transport data for specified segments
    # avoid "efficiency" and "load factor" columns
    car_cols = transport_data.columns[
        ~transport_data.columns.str.contains("efficiency|load factor")
    ]
    nodal_transport_data[car_cols] = nodal_transport_data[car_cols].mul(
        pop_layout["fraction"], axis=0
    )
    # fill missing stats with average data
    stats = [
        "average fuel efficiency",
        "load factor Rail passenger",
        "load factor Rail freight",
        "load factor Heavy goods vehicles",
        "load factor Motor coaches, buses and trolley buses",
    ]
    for stat in stats:
        nodal_transport_data.loc[
            nodal_transport_data[stat] == 0.0,
            stat,
        ] = transport_data[stat].mean()

    return nodal_transport_data


def get_shape(traffic_fn):
    # averaged weekly counts from the year 2010-2015
    traffic = pd.read_csv(traffic_fn, skiprows=2, usecols=["count"]).squeeze("columns")

    # create annual profile take account time zone + summer time
    transport_shape = generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=traffic.values,
    )
    transport_shape = transport_shape / transport_shape.sum()

    return transport_shape


def build_transport_demand(
    traffic_fn_pc,
    traffic_fn_ptw,
    traffic_fn_bus,
    traffic_fn_lcv,
    traffic_fn_hgv,
    airtemp_fn,
    nodes,
    nodal_transport_data,
):
    """
    Returns transport demand per bus in unit km driven [100 km] for the five distinguishable motor vehicle types:
    - "pkw" - pc = passenger cars,
    - "mot" - ptw = powered two-wheelers,
    - "bus" = buses/coaches,
    - "lcv" = light commercial vehicles,
    - "hgv" = heavy goods vehicles.
    """
    # get transport shape per vehicle type
    transport_shape_pc = get_shape(traffic_fn_pc)
    transport_shape_ptw = get_shape(traffic_fn_ptw)
    transport_shape_bus = get_shape(traffic_fn_bus)
    transport_shape_lcv = get_shape(traffic_fn_lcv)
    transport_shape_hgv = get_shape(traffic_fn_hgv)

    # non-electrified rail share
    non_elec_rail = 1 - (
        pop_weighted_energy_totals["electricity rail"]
        / pop_weighted_energy_totals["total rail"]
    )

    # total demand of driven vehicle-km [mio km]
    pc = nodal_transport_data["mio km-driven Passenger cars"]
    ptw = nodal_transport_data["mio km-driven Powered two-wheelers"]
    lcv = nodal_transport_data["mio km-driven Light commercial vehicles"]
    bus = nodal_transport_data[
        "mio km-driven Motor coaches, buses and trolley buses"
    ] + non_elec_rail * nodal_transport_data["mio km-driven Rail passenger"] * (
        nodal_transport_data["load factor Rail passenger"]
        / nodal_transport_data["load factor Motor coaches, buses and trolley buses"]
    )
    hgv = nodal_transport_data[
        "mio km-driven Heavy goods vehicles"
    ] + non_elec_rail * nodal_transport_data["mio km-driven Rail freight"] * (
        nodal_transport_data["load factor Rail freight"]
        / nodal_transport_data["load factor Heavy goods vehicles"]
    )

    def get_demand(profile, total, nyears, seg):
        """
        Returns from total demand [mio km] and given profile
        demand time-series in unit [100 km].
        """

        demand = profile.multiply(total) * 1e4 * nyears

        return pd.concat([demand], keys=[seg], axis=1)

    demand_pc = get_demand(transport_shape_pc, pc, nyears, seg="pc")
    demand_ptw = get_demand(transport_shape_ptw, ptw, nyears, seg="ptw")
    demand_bus = get_demand(transport_shape_bus, bus, nyears, seg="bus")
    demand_lcv = get_demand(transport_shape_lcv, lcv, nyears, seg="lcv")
    demand_hgv = get_demand(transport_shape_hgv, hgv, nyears, seg="hgv")

    return pd.concat(
        [demand_pc, demand_ptw, demand_lcv, demand_hgv, demand_bus], axis=1
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


def bev_availability_profile(
    fn_pc, fn_ptw, fn_lcv, fn_hgv, fn_bus, snapshots, nodes, options
):
    """
    Derive plugged-in availability for electric vehicles.
    """
    # car count in typical week
    traffic_pc = pd.read_csv(fn_pc, skiprows=2, usecols=["count"]).squeeze("columns")
    traffic_ptw = pd.read_csv(fn_ptw, skiprows=2, usecols=["count"]).squeeze("columns")
    traffic_bus = pd.read_csv(fn_bus, skiprows=2, usecols=["count"]).squeeze("columns")
    traffic_lcv = pd.read_csv(fn_lcv, skiprows=2, usecols=["count"]).squeeze("columns")
    traffic_hgv = pd.read_csv(fn_hgv, skiprows=2, usecols=["count"]).squeeze("columns")

    def get_avail(traffic, seg):
        # maximum share plugged-in availability for respective segment
        avail_max = options["bev_avail_max"][seg]
        # average share plugged-in availability for respective segment
        avail_mean = options["bev_avail_mean"][seg]

        # linear scaling, highest when traffic is lowest, decreases if traffic increases
        avail = avail_max - (avail_max - avail_mean) * (traffic - traffic.min()) / (
            traffic.mean() - traffic.min()
        )

        if not avail[avail < 0].empty:
            logger.warning(
                "The BEV availability weekly profile has negative values which can "
                "lead to infeasibility."
            )

        avail_periodic = generate_periodic_profiles(
            dt_index=snapshots,
            nodes=nodes,
            weekly_profile=avail.values,
        )

        return pd.concat([avail_periodic], keys=[seg], axis=1)

    return pd.concat(
        [
            get_avail(traffic_pc, seg="pc"),
            get_avail(traffic_ptw, seg="ptw"),
            get_avail(traffic_bus, seg="bus"),
            get_avail(traffic_lcv, seg="lcv"),
            get_avail(traffic_hgv, seg="hgv"),
        ],
        axis=1,
    )


def bev_dsm_profile(snapshots, nodes, options):

    def get_dsm(seg):
        dsm_week = np.zeros((24 * 7,))

        # assuming that at a certain time ("bev_dsm_restriction_time") EVs have to
        # be charged to a minimum value (defined in bev_dsm_restriction_value)
        dsm_week[
            (np.arange(0, 7, 1) * 24 + options["bev_dsm_restriction_time"][seg])
        ] = options["bev_dsm_restriction_value"][seg]

        dsm_periodic = generate_periodic_profiles(
            dt_index=snapshots,
            nodes=nodes,
            weekly_profile=dsm_week,
        )

        return pd.concat([dsm_periodic], keys=[seg], axis=1)

    return pd.concat(
        [
            get_dsm(seg="pc"),
            get_dsm(seg="ptw"),
            get_dsm(seg="bus"),
            get_dsm(seg="lcv"),
            get_dsm(seg="hgv"),
        ],
        axis=1,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_transport_demand", clusters=128)
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

    n = pypsa.Network(snakemake.input.network)
    nyears = n.snapshot_weightings.generators.sum() / 8760.0

    energy_totals_year = snakemake.params.energy_totals_year
    nodal_transport_data = build_nodal_transport_data(
        snakemake.input.transport_data, pop_layout, energy_totals_year
    )

    transport_demand = build_transport_demand(
        snakemake.input.traffic_data_Pkw,
        snakemake.input.traffic_data_Mot,
        snakemake.input.traffic_data_Bus,
        snakemake.input.traffic_data_Lfw,
        snakemake.input.traffic_data_Lkw,
        snakemake.input.temp_air_total,
        nodes,
        nodal_transport_data,
    )

    avail_profile = bev_availability_profile(
        snakemake.input.traffic_data_Pkw,
        snakemake.input.traffic_data_Mot,
        snakemake.input.traffic_data_Bus,
        snakemake.input.traffic_data_Lfw,
        snakemake.input.traffic_data_Lkw,
        snapshots,
        nodes,
        options,
    )

    dsm_profile = bev_dsm_profile(snapshots, nodes, options)

    nodal_transport_data.to_csv(snakemake.output.transport_data)
    transport_demand.to_csv(snakemake.output.transport_demand)
    avail_profile.to_csv(snakemake.output.avail_profile)
    dsm_profile.to_csv(snakemake.output.dsm_profile)
