# -*- coding: utf-8 -*-
"""
Solve operations network.
"""


import logging

import numpy as np
import pypsa
from helper import override_component_attrs
from solve_network import prepare_network, solve_network

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def set_parameters_from_optimized(n, n_optim):
    lines_typed_i = n.lines.index[n.lines.type != ""]
    n.lines.loc[lines_typed_i, "num_parallel"] = n_optim.lines["num_parallel"].reindex(
        lines_typed_i, fill_value=0.0
    )
    n.lines.loc[lines_typed_i, "s_nom"] = (
        np.sqrt(3)
        * n.lines["type"].map(n.line_types.i_nom)
        * n.lines.bus0.map(n.buses.v_nom)
        * n.lines.num_parallel
    )

    lines_untyped_i = n.lines.index[n.lines.type == ""]
    for attr in ("s_nom", "r", "x"):
        n.lines.loc[lines_untyped_i, attr] = n_optim.lines[attr].reindex(
            lines_untyped_i, fill_value=0.0
        )
    n.lines["s_nom_extendable"] = False

    links_dc_i = n.links.index[n.links.p_nom_extendable]
    n.links.loc[links_dc_i, "p_nom"] = n_optim.links["p_nom_opt"].reindex(
        links_dc_i, fill_value=0.0
    )
    n.links.loc[links_dc_i, "p_nom_extendable"] = False

    gen_extend_i = n.generators.index[n.generators.p_nom_extendable]
    n.generators.loc[gen_extend_i, "p_nom"] = n_optim.generators["p_nom_opt"].reindex(
        gen_extend_i, fill_value=0.0
    )
    n.generators.loc[gen_extend_i, "p_nom_extendable"] = False

    stor_units_extend_i = n.storage_units.index[n.storage_units.p_nom_extendable]
    n.storage_units.loc[stor_units_extend_i, "p_nom"] = n_optim.storage_units[
        "p_nom_opt"
    ].reindex(stor_units_extend_i, fill_value=0.0)
    n.storage_units.loc[stor_units_extend_i, "p_nom_extendable"] = False

    stor_extend_i = n.stores.index[n.stores.e_nom_extendable]
    n.stores.loc[stor_extend_i, "e_nom"] = n_optim.stores["e_nom_opt"].reindex(
        stor_extend_i, fill_value=0.0
    )
    n.stores.loc[stor_extend_i, "e_nom_extendable"] = False

    return n


def remove_unused_components(n, threshold=50):
    logger.info("Remove assets that are barely used to speed things up.")

    for c in n.iterate_components({"Store", "Link", "Generator"}):
        attr = "e_nom" if c.name == "Store" else "p_nom"
        to_remove = c.df.loc[c.df[attr] < threshold].index
        logger.info(f"Removing barely used {c.name}s:\n{to_remove}")
        n.mremove(c.name, to_remove)

    return n


def add_load_shedding(n, voll=1e4):
    logger.info("Add load shedding to all buses.")

    if "load" in n.generators.carrier.unique():
        to_remove = n.generators.query("carrier == 'load'").index
        logger.info(f"Removing pre-existing load shedding:\n{to_remove}")
        n.mremove("Generator", to_remove)

    n.madd(
        "Generator",
        n.buses.index,
        suffix=" load",
        bus=n.buses.index,
        carrier="load",
        marginal_cost=voll,
        p_nom=1e6,
    )

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helper import mock_snakemake

        snakemake = mock_snakemake(
            "solve_operations_network",
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

    overrides = override_component_attrs(snakemake.input.overrides)

    n = pypsa.Network(snakemake.input.pre, override_component_attrs=overrides)
    n_post = pypsa.Network(snakemake.input.post, override_component_attrs=overrides)
    n = set_parameters_from_optimized(n, n_post)
    del n_post

    n = remove_unused_components(n)
    n = add_load_shedding(n)

    opts = snakemake.wildcards.sector_opts.split("-")
    solve_opts = snakemake.config["solving"]["options"]
    solve_opts["skip_iterations"] = True

    n = prepare_network(n, solve_opts)

    n = solve_network(
        n,
        config=snakemake.config,
        opts=opts,
        solver_dir=tmpdir,
        solver_logfile=snakemake.log.solver,
    )

    n.export_to_netcdf(snakemake.output[0])
