# -*- coding: utf-8 -*-
"""
Solve myopic operations network.
"""


import logging

import pandas as pd
import pypsa
from helper import override_component_attrs
from solve_network import prepare_network, solve_network
from solve_operations_network import (
    add_load_shedding,
    remove_unused_components,
    set_parameters_from_optimized,
)

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def prepare_myopic(n, config, store_soc, storage_unit_soc):
    n.stores.e_cyclic = False
    n.storage_units.cyclic_state_of_charge = False

    biomass_stores = n.stores.carrier.str.isin(["solid biomass", "biogas"])
    biomass_potential = n.stores.loc[biomass_stores, "e_initial"]

    # storage level contiguity across years
    n.stores.e_initial = store_soc
    n.storage_units.state_of_charge_initial = storage_unit_soc

    # replace co2 limit with co2 price
    n.remove("GlobalConstraint", "CO2Limit")
    n.stores.at["co2 atmosphere", "marginal_cost"] = -config["co2_price"]

    # handle co2 sequestration
    assert (
        sum(n.stores.carriers == "co2 stored") == 1
    ), "Myopic operation not implemented for spatially resolved CO2 sequestration."
    n.stores.at["co2 stored", "e_nom"] = config["co2_sequestration_limit"] * 1e6  # t/a

    # reset co2 emissions
    n.stores.loc[n.stores.carrier == "co2 stored", "e_initial"] = 0.0
    n.stores.at["co2 atmosphere", "e_initial"] = 0.0

    # replenish fossil gas and oil with 1000 TWh each
    fossil_stores = n.stores.carrier.str.isin(["gas", "oil"])
    n.stores.loc[fossil_stores, "e_initial"] = 1e9
    n.stores.loc[fossil_stores, "e_nom"] = 10e9

    # replenish annual solid biomass and biogas potentials
    n.stores.loc[biomass_stores, "e_initial"] = biomass_potential

    # set storage bidding prices
    bidding_prices = config["bidding_prices"]
    for c in n.iterate_components({"Store", "Link", "StorageUnit"}):
        c.df.marginal_cost.update(c.df.carrier.map(bidding_prices).dropna())

    # deduct industry solid biomass
    assert (
        sum(n.stores.carriers == "solid biomass") == 1
    ), "Myopic operation not implemented for spatially resolved solid biomass."
    n.stores.at["EU solid biomass", "e_initial"] -= (
        n.loads.at["solid biomass for industry", "p_set"] * 8760
    )
    n.remove("Load", "solid biomass for industry")

    return n


def solve_network_myopic(n, config, opts="", **kwargs):
    rolling_horizon = config["operations"]["rolling_horizon"]

    freq = int(pd.infer_freq(n.snapshots)[:-1])
    window = rolling_horizon["window"] * 24 // freq
    overlap = rolling_horizon["overlap"] * 24 // freq
    kept = window - overlap
    length = len(n.snapshots)

    assert (
        kept > 0
    ), f"Overlap ({overlap} days) must be smaller than windows ({window} days)."

    for i in range(length // kept):
        snapshots = n.snapshots[i * kept : (i + 1) * kept + overlap]
        logger.info(f"Optimising operations from {snapshots[0]} to {snapshots[-1]}")

        n = solve_network(n, config, opts=opts, snapshots=snapshots, **kwargs)

        last_kept = n.snapshots[(i + 1) * kept - 1]
        logger.info(f"Setting initial SOCs from {last_kept} for next iteration.\n")

        n.stores.e_initial = n.stores_t.e.loc[last_kept]
        n.storage_units.state_of_charge_initial = n.storage_units_t.state_of_charge.loc[
            last_kept
        ]

    # final segment until end of year
    snapshots = n.snapshots[(i + 1) * kept :]
    n = solve_network(n, config, opts=opts, snapshots=snapshots, **kwargs)

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helper import mock_snakemake

        snakemake = mock_snakemake(
            "solve_operations_network_myopic",
            capacity_year=1952,
            simpl="",
            opts="",
            clusters=37,
            lv=2.0,
            sector_opts="Co2L0-25H-T-H-B-I-A",
            planning_horizons=2030,
            weather_year=2013,
        )

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    tmpdir = snakemake.config["solving"].get("tmpdir")
    if tmpdir is not None:
        from pathlib import Path

        Path(tmpdir).mkdir(parents=True, exist_ok=True)

    config = snakemake.config["operations"]
    overrides = override_component_attrs(snakemake.input.overrides)

    n = pypsa.Network(snakemake.input.pre, override_component_attrs=overrides)

    n_post = pypsa.Network(snakemake.input.post, override_component_attrs=overrides)
    n = set_parameters_from_optimized(n, n_post)
    del n_post

    n_previous = pypsa.Network(
        snakemake.input.previous, override_component_attrs=overrides
    )
    store_soc = n_previous.stores_t.e.iloc[-1]
    storage_unit_soc = n_previous.storage_units_t.state_of_charge.iloc[-1]
    del n_previous

    n = remove_unused_components(n)
    n = add_load_shedding(n)
    n = prepare_myopic(n, config, store_soc, storage_unit_soc)

    opts = snakemake.wildcards.sector_opts.split("-")
    solve_opts = snakemake.config["solving"]["options"]
    solve_opts["skip_iterations"] = True

    n = prepare_network(n, solve_opts)

    n = solve_network_myopic(
        n,
        config=snakemake.config,
        opts=opts,
        solver_dir=tmpdir,
        solver_logfile=snakemake.log.solver,
    )

    n.export_to_netcdf(snakemake.output[0])
