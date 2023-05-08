# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Solves optimal operation and capacity for a network with the option to
iteratively optimize while updating line reactances.

This script is used for optimizing the electrical network as well as the
sector coupled network.

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.

The optimization is based on the :func:`network.optimize` function.
Additionally, some extra constraints specified in :mod:`solve_network` are added.

.. note::

    The rules ``solve_elec_networks`` and ``solve_sector_networks`` run
    the workflow for all scenarios in the configuration file (``scenario:``)
    based on the rule :mod:`solve_network`.
"""
import logging
import re

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    override_component_attrs,
    update_config_with_sector_opts,
)
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)
from pypsa.descriptors import get_switchable_as_dense as get_as_dense


def add_land_use_constraint(n, config):
    if "m" in snakemake.wildcards.clusters:
        _add_land_use_constraint_m(n, config)
    else:
        _add_land_use_constraint(n, config)


def _add_land_use_constraint(n, config):
    # warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in ["solar", "onwind", "offwind-ac", "offwind-dc"]:
        ext_i = (n.generators.carrier == carrier) & ~n.generators.p_nom_extendable
        existing = (
            n.generators.loc[ext_i, "p_nom"]
            .groupby(n.generators.bus.map(n.buses.location))
            .sum()
        )
        existing.index += " " + carrier + "-" + snakemake.wildcards.planning_horizons
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    # check if existing capacities are larger than technical potential
    existing_large = n.generators[
        n.generators["p_nom_min"] > n.generators["p_nom_max"]
    ].index
    if len(existing_large):
        logger.warning(
            f"Existing capacities larger than technical potential for {existing_large},\
                       adjust technical potential to existing capacities"
        )
        n.generators.loc[existing_large, "p_nom_max"] = n.generators.loc[
            existing_large, "p_nom_min"
        ]

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def _add_land_use_constraint_m(n, config):
    # if generators clustering is lower than network clustering, land_use accounting is at generators clusters

    planning_horizons = config["scenario"]["planning_horizons"]
    grouping_years = config["existing_capacities"]["grouping_years"]
    current_horizon = snakemake.wildcards.planning_horizons

    for carrier in ["solar", "onwind", "offwind-ac", "offwind-dc"]:
        existing = n.generators.loc[n.generators.carrier == carrier, "p_nom"]
        ind = list(
            set(
                [
                    i.split(sep=" ")[0] + " " + i.split(sep=" ")[1]
                    for i in existing.index
                ]
            )
        )

        previous_years = [
            str(y)
            for y in planning_horizons + grouping_years
            if y < int(snakemake.wildcards.planning_horizons)
        ]

        for p_year in previous_years:
            ind2 = [
                i for i in ind if i + " " + carrier + "-" + p_year in existing.index
            ]
            sel_current = [i + " " + carrier + "-" + current_horizon for i in ind2]
            sel_p_year = [i + " " + carrier + "-" + p_year for i in ind2]
            n.generators.loc[sel_current, "p_nom_max"] -= existing.loc[
                sel_p_year
            ].rename(lambda x: x[:-4] + current_horizon)

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def add_co2_sequestration_limit(n, limit=200):
    """
    Add a global constraint on the amount of Mt CO2 that can be sequestered.
    """
    n.carriers.loc["co2 stored", "co2_absorptions"] = -1
    n.carriers.co2_absorptions = n.carriers.co2_absorptions.fillna(0)

    limit = limit * 1e6
    for o in opts:
        if "seq" not in o:
            continue
        limit = float(o[o.find("seq") + 3 :]) * 1e6
        break

    n.add(
        "GlobalConstraint",
        "co2_sequestration_limit",
        sense="<=",
        constant=limit,
        type="primary_energy",
        carrier_attribute="co2_absorptions",
    )


def prepare_network(n, solve_opts=None, config=None):
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,  # TODO: check if this can be removed
            n.storage_units_t.inflow,
        ):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    load_shedding = solve_opts.get("load_shedding")
    if load_shedding:
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        # TODO: retrieve color and nice name from config
        n.add("Carrier", "load", color="#dd2e23", nice_name="Load shedding")
        buses_i = n.buses.query("carrier == 'AC'").index
        if not np.isscalar(load_shedding):
            # TODO: do not scale via sign attribute (use Eur/MWh instead of Eur/kWh)
            load_shedding = 1e2  # Eur/kWh

        n.madd(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,  # Eur/kWh
            p_nom=1e9,  # kW
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if config["foresight"] == "myopic":
        add_land_use_constraint(n, config)

    if n.stores.carrier.eq("co2 stored").any():
        limit = config["sector"].get("co2_sequestration_potential", 200)
        add_co2_sequestration_limit(n, limit=limit)

    return n


def add_CCL_constraints(n, config):
    """
    Add CCL (country & carrier limit) constraint to the network.

    Add minimum and maximum levels of generator nominal capacity per carrier
    for individual countries. Opts and path for agg_p_nom_minmax.csv must be defined
    in config.yaml. Default file is available at data/agg_p_nom_minmax.csv.

    Parameters
    ----------
    n : pypsa.Network
    config : dict

    Example
    -------
    scenario:
        opts: [Co2L-CCL-24H]
    electricity:
        agg_p_nom_limits: data/agg_p_nom_minmax.csv
    """
    agg_p_nom_minmax = pd.read_csv(
        config["electricity"]["agg_p_nom_limits"], index_col=[0, 1]
    )
    logger.info("Adding generation capacity constraints per carrier and country")
    p_nom = n.model["Generator-p_nom"]

    gens = n.generators.query("p_nom_extendable").rename_axis(index="Generator-ext")
    grouper = [gens.bus.map(n.buses.country), gens.carrier]
    grouper = xr.DataArray(pd.MultiIndex.from_arrays(grouper), dims=["Generator-ext"])
    lhs = p_nom.groupby(grouper).sum().rename(bus="country")

    minimum = xr.DataArray(agg_p_nom_minmax["min"].dropna()).rename(dim_0="group")
    index = minimum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) >= minimum.loc[index], name="agg_p_nom_min"
        )

    maximum = xr.DataArray(agg_p_nom_minmax["max"].dropna()).rename(dim_0="group")
    index = maximum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) <= maximum.loc[index], name="agg_p_nom_max"
        )


def add_EQ_constraints(n, level, by_country, config, scaling=1e-1):
    """
    Add equity constraints to the network.

    The equity option specifies a certain level x of equity (as in
    EQx) where x is a number between 0 and 1. When this option is set,
    each node in the network is required to produce at least a
    fraction of x of its energy demand locally. For example, when
    EQ0.7 is set, each node is required to produce at least 70% of its
    energy demand locally.

    Locally produced energy includes local renewable generation and
    (some) conventional generation (nuclear, coal, geothermal).
    How conventional generation is dealt with depends on whether the
    model is run in electricity-only mode or is sector-coupled. In
    electricity-only mode, all conventional generation is considered
    local production.

    In the sector-coupled model, however, gas and oil are endogenously
    modelled. Oil is not spatially resolved, meaning that any use of
    oil is considered an import, but any production of oil
    (Fischer-Tropsch) is considered an export. When gas is not
    spatially resolved, it functions the same. When, however, a gas
    network is modelled, imports and exports are instead calculated
    using gas network flow. For now, even locally extracted natural
    gas is considered "imported" for the purposes of this equation.

    For other conventional generation (coal, oil, nuclear) the fuel is
    not endogenously modelled, and this generation is considered local
    (even though implementation-wise nodes have to "import" the fuel
    from a copper-plated "EU" node).

    Optionally the EQ option may be suffixed by the letter "c", which
    makes the equity constraint act on a country level instead of a
    node level.

    Regardless, the equity constraint is only enforced on average over
    the whole year.

    In a sector-coupled network, energy production is generally
    greater than consumption because of efficiency losses in energy
    conversions such as hydrogen production (whereas heat pumps
    actually have an "efficiency" greater than 1). Ignoring these
    losses would lead to a weakening of the equity constraint (i.e. if
    1.5TWh of electricity needs to be produced to satisfy a final
    demand of 1 TWh of energy, even an equity constraint of 100% would
    be satisfied if 1TWh of electricity is produced locally).
    Therefore, for the purpose of the equity constraint, efficiency
    losses in a sector-coupled network are effectively counted as
    demand, and the equity constraint is enforced on the sum of final
    demand and efficiency losses.

    Again in the sector-coupled model, some energy supply and energy
    demand are copper-plated, meaning that they are not spatially
    resolved by only modelled europe-wide. For the purpose of the
    equity constraint in a sector-coupled model, energy supplied from
    a copper-plated carrier (supplied from the "european node") is
    counted as imported, not locally produced. Similarly, energy
    demand for a copper-plated carrier (demanded at the "european
    node") is counted as exported, not locally demanded.

    Parameters
    ----------
    n : pypsa.Network
    o : str

    Example
    -------
    scenario:
        opts: [Co2L-EQ0.7-24H]
    """
    # TODO: Does this work for myopic optimisation?

    # Implementation note: while the equity constraint nominally
    # specifies that a minimum fraction demand be produced locally, in
    # the implementation we enforce a minimum ratio between local
    # production and net energy imports. This because
    #     local_production + net_imports = demand + efficiency_losses
    # for each node (or country), so we can convert a constraint of the form
    #     local_production >= level * (demand + efficiency_losses)
    # to the equivalent:
    #     net_imports <= (1 / level - 1) * local_production
    # or, putting all variables on the right hand side and constants
    # on the left hand side:
    #     (1 - 1 / level) * local_production + net_imports <= 0
    #
    # While leading to an equivalent constraint, we choose this
    # implementation because it allows us to avoid having to calculate
    # efficiency losses explicitly; local production and net imports
    # are slightly easier to deal with.
    #
    # Notes on specific technologies. Locally produced energy comes
    # from the following sources:
    # - Variable renewables (various types of solar, onwind, offwind,
    #   etc.), implemented as generators.
    # - Conventional sources (gas, coal, nuclear), implemented as
    #   either generators or links (depending on whether or not the
    #   model is sector-coupled).
    # - Biomass, biogass, if spatially resolved, implemented as stores.
    # - Hydro, implemented as storageunits.
    # - Ambient heat used in heat pumps, implemented as links.
    # Imports can come in the following forms:
    # - Electricity (AC & DC), implemented as lines and links.
    # - Gas pipelines, implemented as links.
    # - H2 pipelines, implemented as links.
    # - Gas imports (pipeline, LNG, production), implemented as generators.

    if config["foresight"] != "overnight":
        logging.warning(
            "Careful! Equity constraint is only tested for 'overnight' "
            f"foresight models, not '{config['foresight']}' foresight"
        )

    # While we need to group components by bus location in the
    # sector-coupled model, there is no "location" column in the
    # electricity-only model.
    location = (
        n.buses.location
        if "location" in n.buses.columns
        else pd.Series(n.buses.index, index=n.buses.index)
    )

    def group(df, b="bus"):
        """
        Group given dataframe by bus location or country.

        The optional argument `b` allows clustering by bus0 or bus1 for
        lines and links.
        """
        if by_country:
            return df[b].map(location).map(n.buses.country).to_xarray()
        else:
            return df[b].map(location).to_xarray()

    # Local production by generators. Note: the network may not
    # actually have all these generators (for instance some
    # conventional generators are implemented as links in the
    # sector-coupled model; heating sector might not be turned on),
    # but we list all that might be in the network.
    local_gen_carriers = list(
        set(
            config["electricity"]["extendable_carriers"]["Generator"]
            + config["electricity"]["conventional_carriers"]
            + config["electricity"]["renewable_carriers"]
            + [c for c in n.generators.carrier if "solar thermal" in c]
            + ["solar rooftop", "wave"]
        )
    )
    local_gen_i = n.generators.loc[
        n.generators.carrier.isin(local_gen_carriers)
        & (n.generators.bus.map(location) != "EU")
    ].index
    local_gen_p = (
        n.model["Generator-p"]
        .loc[:, local_gen_i]
        .groupby(group(n.generators.loc[local_gen_i]))
        .sum()
    )
    local_gen = (local_gen_p * n.snapshot_weightings.generators).sum("snapshot")

    # Hydro production; the only local production from a StorageUnit.
    local_hydro_i = n.storage_units.loc[n.storage_units.carrier == "hydro"].index
    local_hydro_p = (
        n.model["StorageUnit-p_dispatch"]
        .loc[:, local_hydro_i]
        .groupby(group(n.storage_units.loc[local_hydro_i]))
        .sum()
    )
    local_hydro = (local_hydro_p * n.snapshot_weightings.stores).sum("snapshot")

    # Biomass and biogas; these are only considered locally produced
    # if spatially resolved, not if they belong to an "EU" node. They
    # are modelled as stores with initial capacity to model a finite
    # yearly supply; the difference between initial and final capacity
    # is the total local production.
    local_bio_i = n.stores.loc[
        n.stores.carrier.isin(["biogas", "solid biomass"])
        & (n.stores.bus.map(location) != "EU")
    ].index
    # Building the following linear expression only works if it's non-empty
    if len(local_bio_i) > 0:
        local_bio_first_e = n.model["Store-e"].loc[n.snapshots[0], local_bio_i]
        local_bio_last_e = n.model["Store-e"].loc[n.snapshots[-1], local_bio_i]
        local_bio_p = local_bio_first_e - local_bio_last_e
        local_bio = local_bio_p.groupby(group(n.stores.loc[local_bio_i])).sum()
    else:
        local_bio = None

    # Conventional generation in the sector-coupled model. These are
    # modelled as links in order to take the CO2 cycle into account.
    # All of these are counted as local production even if the links
    # may take their fuel from an "EU" node, except for gas and oil,
    # which are modelled endogenously and is counted under imports /
    # exports.
    conv_carriers = config["sector"].get("conventional_generation", {})
    conv_carriers = [
        gen for gen, carrier in conv_carriers.items() if carrier not in ["gas", "oil"]
    ]
    if config["sector"].get("coal_cc") and not "coal" in conv_carriers:
        conv_carriers.append("coal")
    local_conv_gen_i = n.links.loc[n.links.carrier.isin(conv_carriers)].index
    if len(local_conv_gen_i) > 0:
        local_conv_gen_p = n.model["Link-p"].loc[:, local_conv_gen_i]
        # These links have efficiencies, which we divide by since we only
        # want to count the _output_ of each conventional generator as
        # local generation for the equity balance.
        efficiencies = n.links.loc[local_conv_gen_i, "efficiency"]
        local_conv_gen_p = (
            (local_conv_gen_p / efficiencies)
            .groupby(group(n.links.loc[local_conv_gen_i], b="bus1"))
            .sum()
            .rename({"bus1": "bus"})
        )
        local_conv_gen = (local_conv_gen_p * n.snapshot_weightings.generators).sum(
            "snapshot"
        )
    else:
        local_conv_gen = None

    # TODO: should we (in prepare_sector_network.py) model gas
    # pipeline imports from outside the EU and LNG imports separately
    # from gas extraction / production? Then we could model gas
    # extraction as locally produced energy.

    # Ambient heat for heat pumps
    heat_pump_i = n.links.filter(like="heat pump", axis="rows").index
    if len(heat_pump_i) > 0:
        # To get the ambient heat extracted, we subtract 1 from the
        # efficiency of the heat pump (where "efficiency" is really COP
        # for heat pumps).
        from_ambient = n.links_t["efficiency"].loc[:, heat_pump_i] - 1
        local_heat_from_ambient_p = n.model["Link-p"].loc[:, heat_pump_i]
        local_heat_from_ambient = (
            (local_heat_from_ambient_p * from_ambient)
            .groupby(group(n.links.loc[heat_pump_i], b="bus1"))
            .sum()
            .rename({"bus1": "bus"})
        )
        local_heat_from_ambient = (
            local_heat_from_ambient * n.snapshot_weightings.generators
        ).sum("snapshot")
    else:
        local_heat_from_ambient = None

    # Total locally produced energy
    local_energy = sum(
        e
        for e in [
            local_gen,
            local_hydro,
            local_bio,
            local_conv_gen,
            local_heat_from_ambient,
        ]
        if e is not None
    )

    # Now it's time to collect imports: electricity, hydrogen & gas
    # pipeline, other gas, biomass, gas terminals & production.

    # Start with net electricity imports.
    lines_cross_region_i = n.lines.loc[
        (group(n.lines, b="bus0") != group(n.lines, b="bus1")).to_numpy()
    ].index
    # Build linear expression representing net imports (i.e. imports -
    # exports) for each bus/country.
    lines_in_s = (
        n.model["Line-s"]
        .loc[:, lines_cross_region_i]
        .groupby(group(n.lines.loc[lines_cross_region_i], b="bus1"))
        .sum()
        .rename({"bus1": "bus"})
    ) - (
        n.model["Line-s"]
        .loc[:, lines_cross_region_i]
        .groupby(group(n.lines.loc[lines_cross_region_i], b="bus0"))
        .sum()
        .rename({"bus0": "bus"})
    )
    line_imports = (lines_in_s * n.snapshot_weightings.generators).sum("snapshot")

    # Link net imports, representing all net energy imports of various
    # carriers that are implemented as links. We list all possible
    # link carriers that could be represented in the network; some
    # might not be present in some networks depending on the sector
    # configuration. Note that we do not count efficiencies here (e.g.
    # for oil boilers that import oil) since efficiency losses are
    # counted as "local demand".
    link_import_carriers = [
        # Pipeline imports / exports
        "H2 pipeline",
        "H2 pipeline retrofitted",
        "gas pipeline",
        "gas pipeline new",
        # Solid biomass
        "solid biomass transport",
        # DC electricity
        "DC",
        # Oil (imports / exports between spatial nodes and "EU" node)
        "Fischer-Tropsch",
        "biomass to liquid",
        "residential rural oil boiler",
        "services rural oil boiler",
        "residential urban decentral oil boiler",
        "services urban decentral oil boiler",
        "oil",  # Oil powerplant (from `prepare_sector_network.add_generation`)
        # Gas (imports / exports between spatial nodes and "EU" node,
        # only cross-region if gas is not spatially resolved)
        "Sabatier",
        "helmeth",
        "SMR CC",
        "SMR",
        "biogas to gas",
        "BioSNG",
        "residential rural gas boiler",
        "services rural gas boiler",
        "residential urban decentral gas boiler",
        "services urban decentral gas boiler",
        "urban central gas boiler",
        "urban central gas CHP",
        "urban central gas CHP CC",
        "residential rural micro gas CHP",
        "services rural micro gas CHP",
        "residential urban decentral micro gas CHP",
        "services urban decentral micro gas CHP",
        "allam",
        "OCGT",
        "CCGT",
    ]
    links_cross_region_i = (
        n.links.loc[(group(n.links, b="bus0") != group(n.links, b="bus1")).to_numpy()]
        .loc[n.links.carrier.isin(link_import_carriers)]
        .index
    )
    # Build linear expression representing net imports (i.e. imports -
    # exports) for each bus/country.
    links_in_p = (
        n.model["Link-p"]
        .loc[:, links_cross_region_i]
        .groupby(group(n.links.loc[links_cross_region_i], b="bus1"))
        .sum()
        .rename({"bus1": "bus"})
    ) - (
        n.model["Link-p"]
        .loc[:, links_cross_region_i]
        .groupby(group(n.links.loc[links_cross_region_i], b="bus0"))
        .sum()
        .rename({"bus0": "bus"})
    )
    link_imports = (links_in_p * n.snapshot_weightings.generators).sum("snapshot")

    # Gas imports by pipeline from outside of Europe, LNG terminal or
    # gas production (all modelled as generators).
    gas_import_i = n.generators.loc[n.generators.carrier == "gas"].index
    if len(gas_import_i) > 0:
        gas_import_p = (
            n.model["Generator-p"]
            .loc[:, gas_import_i]
            .groupby(group(n.generators.loc[gas_import_i]))
            .sum()
        )
        gas_imports = (gas_import_p * n.snapshot_weightings.generators).sum("snapshot")
    else:
        gas_imports = None

    imported_energy = sum(
        i for i in [line_imports, link_imports, gas_imports] if i is not None
    )

    local_factor = 1 - 1 / level

    n.model.add_constraints(
        local_factor * local_energy + imported_energy <= 0, name="equity_min"
    )


def add_BAU_constraints(n, config):
    """
    Add a per-carrier minimal overall capacity.

    BAU_mincapacities and opts must be adjusted in the config.yaml.

    Parameters
    ----------
    n : pypsa.Network
    config : dict

    Example
    -------
    scenario:
        opts: [Co2L-BAU-24H]
    electricity:
        BAU_mincapacities:
            solar: 0
            onwind: 0
            OCGT: 100000
            offwind-ac: 0
            offwind-dc: 0
    Which sets minimum expansion across all nodes e.g. in Europe to 100GW.
    OCGT bus 1 + OCGT bus 2 + ... > 100000
    """
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    p_nom = n.model["Generator-p_nom"]
    ext_i = n.generators.query("p_nom_extendable")
    ext_carrier_i = xr.DataArray(ext_i.carrier.rename_axis("Generator-ext"))
    lhs = p_nom.groupby(ext_carrier_i).sum()
    index = mincaps.index.intersection(lhs.indexes["carrier"])
    rhs = mincaps[index].rename_axis("carrier")
    n.model.add_constraints(lhs >= rhs, name="bau_mincaps")


# TODO: think about removing or make per country
def add_SAFE_constraints(n, config):
    """
    Add a capacity reserve margin of a certain fraction above the peak demand.
    Renewable generators and storage do not contribute. Ignores network.

    Parameters
    ----------
        n : pypsa.Network
        config : dict

    Example
    -------
    config.yaml requires to specify opts:

    scenario:
        opts: [Co2L-SAFE-24H]
    electricity:
        SAFE_reservemargin: 0.1
    Which sets a reserve margin of 10% above the peak demand.
    """
    peakdemand = n.loads_t.p_set.sum(axis=1).max()
    margin = 1.0 + config["electricity"]["SAFE_reservemargin"]
    reserve_margin = peakdemand * margin
    # TODO: do not take this from the plotting config!
    conv_techs = config["plotting"]["conv_techs"]
    ext_gens_i = n.generators.query("carrier in @conv_techs & p_nom_extendable").index
    p_nom = n.model["Generator-p_nom"].loc[ext_gens_i]
    lhs = p_nom.sum()
    exist_conv_caps = n.generators.query(
        "~p_nom_extendable & carrier in @conv_techs"
    ).p_nom.sum()
    rhs = reserve_margin - exist_conv_caps
    n.model.add_constraints(lhs >= rhs, name="safe_mintotalcap")


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.

    Parameters
    ----------
        n : pypsa.Network
        sns: pd.DatetimeIndex
        config : dict

    Example:
    --------
    config.yaml requires to specify operational_reserve:
    operational_reserve: # like https://genxproject.github.io/GenX/dev/core/#Reserves
        activate: true
        epsilon_load: 0.02 # percentage of load at each snapshot
        epsilon_vres: 0.02 # percentage of VRES at each snapshot
        contingency: 400000 # MW
    """
    reserve_config = config["electricity"]["operational_reserve"]
    EPSILON_LOAD = reserve_config["epsilon_load"]
    EPSILON_VRES = reserve_config["epsilon_vres"]
    CONTINGENCY = reserve_config["contingency"]

    # Reserve Variables
    n.model.add_variables(
        0, np.inf, coords=[sns, n.generators.index], name="Generator-r"
    )
    reserve = n.model["Generator-r"]
    summed_reserve = reserve.sum("Generator")

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        p_nom_vres = (
            n.model["Generator-p_nom"]
            .loc[vres_i.intersection(ext_i)]
            .rename({"Generator-ext": "Generator"})
        )
        lhs = summed_reserve + (p_nom_vres * (-EPSILON_VRES * capacity_factor)).sum(
            "Generator"
        )

    # Total demand per t
    demand = get_as_dense(n, "Load", "p_set").sum(axis=1)

    # VRES potential of non extendable generators
    capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
    renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
    potential = (capacity_factor * renewable_capacity).sum(axis=1)

    # Right-hand-side
    rhs = EPSILON_LOAD * demand + EPSILON_VRES * potential + CONTINGENCY

    n.model.add_constraints(lhs >= rhs, name="reserve_margin")

    # additional constraint that capacity is not exceeded
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = n.model["Generator-p"]
    reserve = n.model["Generator-r"]

    capacity_variable = n.model["Generator-p_nom"].rename(
        {"Generator-ext": "Generator"}
    )
    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    lhs = dispatch + reserve - capacity_variable * p_max_pu[ext_i]

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)

    n.model.add_constraints(lhs <= rhs, name="Generator-p-reserve-upper")


def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.index.str.contains("battery discharger")
    charger_bool = n.links.index.str.contains("battery charger")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


def add_chp_constraints(n):
    electric = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("electric")
    )
    heat = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("heat")
    )

    electric_ext = n.links[electric].query("p_nom_extendable").index
    heat_ext = n.links[heat].query("p_nom_extendable").index

    electric_fix = n.links[electric].query("~p_nom_extendable").index
    heat_fix = n.links[heat].query("~p_nom_extendable").index

    p = n.model["Link-p"]  # dimension: [time, link]

    # output ratio between heat and electricity and top_iso_fuel_line for extendable
    if not electric_ext.empty:
        p_nom = n.model["Link-p_nom"]

        lhs = (
            p_nom.loc[electric_ext]
            * (n.links.p_nom_ratio * n.links.efficiency)[electric_ext].values
            - p_nom.loc[heat_ext] * n.links.efficiency[heat_ext].values
        )
        n.model.add_constraints(lhs == 0, name="chplink-fix_p_nom_ratio")

        rename = {"Link-ext": "Link"}
        lhs = (
            p.loc[:, electric_ext]
            + p.loc[:, heat_ext]
            - p_nom.rename(rename).loc[electric_ext]
        )
        n.model.add_constraints(lhs <= 0, name="chplink-top_iso_fuel_line_ext")

    # top_iso_fuel_line for fixed
    if not electric_fix.empty:
        lhs = p.loc[:, electric_fix] + p.loc[:, heat_fix]
        rhs = n.links.p_nom[electric_fix]
        n.model.add_constraints(lhs <= rhs, name="chplink-top_iso_fuel_line_fix")

    # back-pressure
    if not electric.empty:
        lhs = (
            p.loc[:, heat] * (n.links.efficiency[heat] * n.links.c_b[electric].values)
            - p.loc[:, electric] * n.links.efficiency[electric]
        )
        n.model.add_constraints(lhs <= rhs, name="chplink-backpressure")


def add_pipe_retrofit_constraint(n):
    """
    Add constraint for retrofitting existing CH4 pipelines to H2 pipelines.
    """
    gas_pipes_i = n.links.query("carrier == 'gas pipeline' and p_nom_extendable").index
    h2_retrofitted_i = n.links.query(
        "carrier == 'H2 pipeline retrofitted' and p_nom_extendable"
    ).index

    if h2_retrofitted_i.empty or gas_pipes_i.empty:
        return

    p_nom = n.model["Link-p_nom"]

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    lhs = p_nom.loc[gas_pipes_i] + CH4_per_H2 * p_nom.loc[h2_retrofitted_i]
    rhs = n.links.p_nom[gas_pipes_i].rename_axis("Link-ext")

    n.model.add_constraints(lhs == rhs, name="Link-pipe_retrofit")


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments ``opts`` and
    ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    if "BAU" in opts and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if "SAFE" in opts and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if "CCL" in opts and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config)
    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)
    for o in opts:
        if "EQ" in o:
            EQ_regex = "EQ(0\.[0-9]+)(c?)"  # Ex.: EQ0.75c
            m = re.search(EQ_regex, o)
            if m is not None:
                level = float(m.group(1))
                by_country = True if m.group(2) == "c" else False
                add_EQ_constraints(n, level, by_country, config)
            else:
                logging.warning(f"Invalid EQ option: {o}")
    add_battery_constraints(n)
    add_pipe_retrofit_constraint(n)


def solve_network(n, config, opts="", **kwargs):
    set_of_options = config["solving"]["solver"]["options"]
    solver_options = (
        config["solving"]["solver_options"][set_of_options] if set_of_options else {}
    )
    solver_name = config["solving"]["solver"]["name"]
    cf_solving = config["solving"]["options"]
    track_iterations = cf_solving.get("track_iterations", False)
    min_iterations = cf_solving.get("min_iterations", 4)
    max_iterations = cf_solving.get("max_iterations", 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    skip_iterations = cf_solving.get("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    if skip_iterations:
        status, condition = n.optimize(
            solver_name=solver_name,
            extra_functionality=extra_functionality,
            **solver_options,
            **kwargs,
        )
    else:
        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
            solver_name=solver_name,
            track_iterations=track_iterations,
            min_iterations=min_iterations,
            max_iterations=max_iterations,
            extra_functionality=extra_functionality,
            **solver_options,
            **kwargs,
        )

    if status != "ok":
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'"
        )
    if "infeasible" in condition:
        raise RuntimeError("Solving status 'infeasible'")

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_sector_network",
            configfiles="test/config.overnight.yaml",
            simpl="",
            opts="",
            clusters="5",
            ll="v1.5",
            sector_opts="CO2L0-24H-T-H-B-I-A-solar+p3-dist1",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    if "sector_opts" in snakemake.wildcards.keys():
        update_config_with_sector_opts(
            snakemake.config, snakemake.wildcards.sector_opts
        )

    opts = snakemake.wildcards.opts
    if "sector_opts" in snakemake.wildcards.keys():
        opts += "-" + snakemake.wildcards.sector_opts
    opts = [o for o in opts.split("-") if o != ""]
    solve_opts = snakemake.config["solving"]["options"]

    np.random.seed(solve_opts.get("seed", 123))

    fn = getattr(snakemake.log, "memory", None)
    with memory_logger(filename=fn, interval=30.0) as mem:
        if "overrides" in snakemake.input.keys():
            overrides = override_component_attrs(snakemake.input.overrides)
            n = pypsa.Network(
                snakemake.input.network, override_component_attrs=overrides
            )
        else:
            n = pypsa.Network(snakemake.input.network)

        n = prepare_network(n, solve_opts, config=snakemake.config)

        n = solve_network(
            n, config=snakemake.config, opts=opts, log_fn=snakemake.log.solver
        )

        n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
