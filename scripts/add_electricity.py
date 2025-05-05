# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Adds existing electrical generators, hydro-electric plants as well as
greenfield and battery and hydrogen storage to the clustered network.

Description
-----------

The rule :mod:`add_electricity` ties all the different data inputs from the
preceding rules together into a detailed PyPSA network that is stored in
``networks/base_s_{clusters}_elec.nc``. It includes:

- today's transmission topology and transfer capacities (optionally including
  lines which are under construction according to the config settings ``lines:
  under_construction`` and ``links: under_construction``),
- today's thermal and hydro power generation capacities (for the technologies
  listed in the config setting ``electricity: conventional_carriers``), and
- today's load time-series (upsampled in a top-down approach according to
  population and gross domestic product)

It further adds extendable ``generators`` with **zero** capacity for

- photovoltaic, onshore and AC- as well as DC-connected offshore wind
  installations with today's locational, hourly wind and solar capacity factors
  (but **no** current capacities),
- additional open- and combined-cycle gas turbines (if ``OCGT`` and/or ``CCGT``
  is listed in the config setting ``electricity: extendable_carriers``)

Furthermore, it attaches additional extendable components to the clustered
network with **zero** initial capacity:

- ``StorageUnits`` of carrier 'H2' and/or 'battery'. If this option is chosen,
  every bus is given an extendable ``StorageUnit`` of the corresponding carrier.
  The energy and power capacities are linked through a parameter that specifies
  the energy capacity as maximum hours at full dispatch power and is configured
  in ``electricity: max_hours:``. This linkage leads to one investment variable
  per storage unit. The default ``max_hours`` lead to long-term hydrogen and
  short-term battery storage units.

- ``Stores`` of carrier 'H2' and/or 'battery' in combination with ``Links``. If
  this option is chosen, the script adds extra buses with corresponding carrier
  where energy ``Stores`` are attached and which are connected to the
  corresponding power buses via two links, one each for charging and
  discharging. This leads to three investment variables for the energy capacity,
  charging and discharging capacity of the storage unit.
"""

import logging
from collections.abc import Iterable
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from pypsa.clustering.spatial import DEFAULT_ONE_PORT_STRATEGIES, normed_or_uniform

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    rename_techs,
    set_scenario_config,
    update_p_nom_max,
)

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(s: pd.Series) -> pd.Series:
    """
    Normalize a pandas Series by dividing each element by the sum of all elements.

    Parameters
    ----------
    s : pd.Series
        Input series to normalize

    Returns
    -------
    pd.Series
        Normalized series where all elements sum to 1
    """
    return s / s.sum()


def flatten(t: Iterable[Any]) -> str:
    return " ".join(map(str, t))


def calculate_annuity(n: float, r: float | pd.Series) -> float | pd.Series:
    """
    Calculate the annuity factor for an asset with lifetime n years and discount rate r.

    The annuity factor is used to calculate the annual payment required to pay off a loan
    over n years at interest rate r. For example, annuity(20, 0.05) * 20 = 1.6.

    Parameters
    ----------
    n : float
        Lifetime of the asset in years
    r : float | pd.Series
        Discount rate (interest rate). Can be a single float or a pandas Series of rates.

    Returns
    -------
    float | pd.Series
        Annuity factor. Returns a float if r is float, or pd.Series if r is pd.Series.

    Examples
    --------
    >>> calculate_annuity(20, 0.05)
    0.08024258718774728
    """
    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def add_missing_carriers(n, carriers):
    """
    Function to add missing carriers to the network without raising errors.
    """
    missing_carriers = set(carriers) - set(n.carriers.index)
    if len(missing_carriers) > 0:
        n.add("Carrier", missing_carriers)


def sanitize_carriers(n, config):
    """
    Sanitize the carrier information in a PyPSA Network object.

    The function ensures that all unique carrier names are present in the network's
    carriers attribute, and adds nice names and colors for each carrier according
    to the provided configuration dictionary.

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA Network object that represents an electrical power system.
    config : dict
        A dictionary containing configuration information, specifically the
        "plotting" key with "nice_names" and "tech_colors" keys for carriers.

    Returns
    -------
    None
        The function modifies the 'n' PyPSA Network object in-place, updating the
        carriers attribute with nice names and colors.

    Warnings
    --------
    Raises a warning if any carrier's "tech_colors" are not defined in the config dictionary.
    """

    for c in n.iterate_components():
        if "carrier" in c.df:
            add_missing_carriers(n, c.df.carrier)

    carrier_i = n.carriers.index
    nice_names = (
        pd.Series(config["plotting"]["nice_names"])
        .reindex(carrier_i)
        .fillna(carrier_i.to_series())
    )
    n.carriers["nice_name"] = n.carriers.nice_name.where(
        n.carriers.nice_name != "", nice_names
    )

    tech_colors = config["plotting"]["tech_colors"]
    colors = pd.Series(tech_colors).reindex(carrier_i)
    # try to fill missing colors with tech_colors after renaming
    missing_colors_i = colors[colors.isna()].index
    colors[missing_colors_i] = missing_colors_i.map(rename_techs).map(tech_colors)
    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(f"tech_colors for carriers {missing_i} not defined in config.")
    n.carriers["color"] = n.carriers.color.where(n.carriers.color != "", colors)


def sanitize_locations(n):
    if "location" in n.buses.columns:
        n.buses["x"] = n.buses.x.where(n.buses.x != 0, n.buses.location.map(n.buses.x))
        n.buses["y"] = n.buses.y.where(n.buses.y != 0, n.buses.location.map(n.buses.y))
        n.buses["country"] = n.buses.country.where(
            n.buses.country.ne("") & n.buses.country.notnull(),
            n.buses.location.map(n.buses.country),
        )


def add_co2_emissions(n, costs, carriers):
    """
    Add CO2 emissions to the network's carriers attribute.
    """
    suptechs = n.carriers.loc[carriers].index.str.split("-").str[0]
    n.carriers.loc[carriers, "co2_emissions"] = costs.loc[
        suptechs, "CO2 intensity"
    ].values


def load_costs(
    cost_file: str, config: dict, max_hours: dict = None, nyears: float = 1.0
) -> pd.DataFrame:
    """
    Load cost data from CSV and prepare it.

    Parameters
    ----------
    cost_file : str
        Path to the CSV file containing cost data
    config : dict
        Dictionary containing cost-related configuration parameters
    max_hours : dict, optional
        Dictionary specifying maximum hours for storage technologies
    nyears : float, optional
        Number of years for investment, by default 1.0

    Returns
    -------
    costs : pd.DataFrame
        DataFrame containing the processed cost data
    """
    # Copy marginal_cost and capital_cost for backward compatibility
    for key in ("marginal_cost", "capital_cost"):
        if key in config:
            config["overwrites"][key] = config[key]

    # set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("/GW"), "value"] /= 1e3

    costs.unit = costs.unit.str.replace("/kW", "/MW")
    costs.unit = costs.unit.str.replace("/GW", "/MW")

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = costs.value.unstack(level=1).groupby("technology").sum(min_count=1)
    costs = costs.fillna(config["fill_values"])

    # Process overwrites for various attributes
    for attr in ("investment", "lifetime", "FOM", "VOM", "efficiency", "fuel"):
        overwrites = config["overwrites"].get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            logger.info(f"Overwriting {attr} with:\n{overwrites}")

    annuity_factor = calculate_annuity(costs["lifetime"], costs["discount rate"])
    annuity_factor_fom = annuity_factor + costs["FOM"] / 100.0
    costs["capital_cost"] = annuity_factor_fom * costs["investment"] * nyears

    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    costs.at["OCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]
    costs.at["CCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]

    costs.at["solar", "capital_cost"] = costs.at["solar-utility", "capital_cost"]
    costs = costs.rename({"solar-utility single-axis tracking": "solar-hsat"})

    # Calculate storage costs if max_hours is provided
    if max_hours is not None:

        def costs_for_storage(store, link1, link2=None, max_hours=1.0):
            capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
            if link2 is not None:
                capital_cost += link2["capital_cost"]
            return pd.Series(
                {
                    "capital_cost": capital_cost,
                    "marginal_cost": 0.0,
                    "CO2 intensity": 0.0,
                }
            )

        costs.loc["battery"] = costs_for_storage(
            costs.loc["battery storage"],
            costs.loc["battery inverter"],
            max_hours=max_hours["battery"],
        )
        costs.loc["H2"] = costs_for_storage(
            costs.loc["hydrogen storage underground"],
            costs.loc["fuel cell"],
            costs.loc["electrolysis"],
            max_hours=max_hours["H2"],
        )

    for attr in ("marginal_cost", "capital_cost"):
        overwrites = config["overwrites"].get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            idx = overwrites.index.intersection(costs.index)
            costs.loc[idx, attr] = overwrites.loc[idx]
            logger.info(f"Overwriting {attr} with:\n{overwrites}")

    return costs


def load_and_aggregate_powerplants(
    ppl_fn: str,
    costs: pd.DataFrame,
    consider_efficiency_classes: bool = False,
    aggregation_strategies: dict = None,
    exclude_carriers: list = None,
) -> pd.DataFrame:
    if not aggregation_strategies:
        aggregation_strategies = {}

    if not exclude_carriers:
        exclude_carriers = []

    carrier_dict = {
        "ocgt": "OCGT",
        "ccgt": "CCGT",
        "bioenergy": "biomass",
        "ccgt, thermal": "CCGT",
        "hard coal": "coal",
    }
    tech_dict = {
        "Run-Of-River": "ror",
        "Reservoir": "hydro",
        "Pumped Storage": "PHS",
    }
    ppl = (
        pd.read_csv(ppl_fn, index_col=0, dtype={"bus": "str"})
        .powerplant.to_pypsa_names()
        .rename(columns=str.lower)
        .replace({"carrier": carrier_dict, "technology": tech_dict})
    )

    # Replace carriers "natural gas" and "hydro" with the respective technology;
    # OCGT or CCGT and hydro, PHS, or ror)
    ppl["carrier"] = ppl.carrier.where(
        ~ppl.carrier.isin(["hydro", "natural gas"]), ppl.technology
    )

    cost_columns = [
        "VOM",
        "FOM",
        "efficiency",
        "capital_cost",
        "marginal_cost",
        "fuel",
        "lifetime",
    ]
    ppl = ppl.join(costs[cost_columns], on="carrier", rsuffix="_r")

    ppl["efficiency"] = ppl.efficiency.combine_first(ppl.efficiency_r)
    ppl["lifetime"] = (ppl.dateout - ppl.datein).fillna(np.inf)
    ppl["build_year"] = ppl.datein.fillna(0).astype(int)
    ppl["marginal_cost"] = (
        ppl.carrier.map(costs.VOM) + ppl.carrier.map(costs.fuel) / ppl.efficiency
    )

    strategies = {
        **DEFAULT_ONE_PORT_STRATEGIES,
        **{"country": "first"},
        **aggregation_strategies.get("generators", {}),
    }
    strategies = {k: v for k, v in strategies.items() if k in ppl.columns}

    to_aggregate = ~ppl.carrier.isin(exclude_carriers)
    df = ppl[to_aggregate].copy()

    if consider_efficiency_classes:
        for c in df.carrier.unique():
            df_c = df.query("carrier == @c")
            low = df_c.efficiency.quantile(0.10)
            high = df_c.efficiency.quantile(0.90)
            if low < high:
                labels = ["low", "medium", "high"]
                suffix = pd.cut(
                    df_c.efficiency, bins=[0, low, high, 1], labels=labels
                ).astype(str)
                df.update({"carrier": df_c.carrier + " " + suffix + " efficiency"})

    grouper = ["bus", "carrier"]
    weights = df.groupby(grouper).p_nom.transform(normed_or_uniform)

    for k, v in strategies.items():
        if v == "capacity_weighted_average":
            df[k] = df[k] * weights
            strategies[k] = pd.Series.sum

    aggregated = df.groupby(grouper, as_index=False).agg(strategies)
    aggregated.index = aggregated.bus + " " + aggregated.carrier
    aggregated.build_year = aggregated.build_year.astype(int)

    disaggregated = ppl[~to_aggregate][aggregated.columns].copy()
    disaggregated.index = (
        disaggregated.bus
        + " "
        + disaggregated.carrier
        + " "
        + disaggregated.index.astype(str)
    )

    return pd.concat([aggregated, disaggregated])


def attach_load(
    n: pypsa.Network,
    load_fn: str,
    busmap_fn: str,
    scaling: float = 1.0,
) -> None:
    """
    Attach load data to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the load data to.
    load_fn : str
        Path to the load data file.
    busmap_fn : str
        Path to the busmap file.
    scaling : float, optional
        Scaling factor for the load data, by default 1.0.
    """
    load = (
        xr.open_dataarray(load_fn).to_dataframe().squeeze(axis=1).unstack(level="time")
    )

    # apply clustering busmap
    busmap = pd.read_csv(busmap_fn, dtype=str).set_index("Bus").squeeze()
    load = load.groupby(busmap).sum().T

    logger.info(f"Load data scaled by factor {scaling}.")
    load *= scaling

    n.add("Load", load.columns, bus=load.columns, p_set=load)  # carrier="electricity"


def set_transmission_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    line_length_factor: float = 1.0,
    link_length_factor: float = 1.0,
) -> None:
    """
    Set the transmission costs for lines and links in the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to set the transmission costs for.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    line_length_factor : float, optional
        Factor to scale the line length, by default 1.0.
    link_length_factor : float, optional
        Factor to scale the link length, by default 1.0.
    """
    n.lines["capital_cost"] = (
        n.lines["length"]
        * line_length_factor
        * costs.at["HVAC overhead", "capital_cost"]
    )

    if n.links.empty:
        return

    dc_b = n.links.carrier == "DC"

    # If there are no dc links, then the 'underwater_fraction' column
    # may be missing. Therefore we have to return here.
    if n.links.loc[dc_b].empty:
        return

    costs = (
        n.links.loc[dc_b, "length"]
        * link_length_factor
        * (
            (1.0 - n.links.loc[dc_b, "underwater_fraction"])
            * costs.at["HVDC overhead", "capital_cost"]
            + n.links.loc[dc_b, "underwater_fraction"]
            * costs.at["HVDC submarine", "capital_cost"]
        )
        + costs.at["HVDC inverter pair", "capital_cost"]
    )
    n.links.loc[dc_b, "capital_cost"] = costs


def attach_wind_and_solar(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    profile_filenames: dict,
    carriers: list | set,
    extendable_carriers: list | set,
    line_length_factor: float = 1.0,
    landfall_lengths: dict = None,
) -> None:
    """
    Attach wind and solar generators to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the generators to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    ppl : pd.DataFrame
        DataFrame containing the power plant data.
    profile_filenames : dict
        Dictionary containing the paths to the wind and solar profiles.
    carriers : list | set
        List of renewable energy carriers to attach.
    extendable_carriers : list | set
        List of extendable renewable energy carriers.
    line_length_factor : float, optional
        Factor to scale the line length, by default 1.0.
    landfall_lengths : dict, optional
        Dictionary containing the landfall lengths for offshore wind, by default None.
    """
    add_missing_carriers(n, carriers)

    if landfall_lengths is None:
        landfall_lengths = {}

    for car in carriers:
        if car == "hydro":
            continue

        landfall_length = landfall_lengths.get(car, 0.0)

        with xr.open_dataset(profile_filenames["profile_" + car]) as ds:
            if ds.indexes["bus"].empty:
                continue

            # if-statement for compatibility with old profiles
            if "year" in ds.indexes:
                ds = ds.sel(year=ds.year.min(), drop=True)

            ds = ds.stack(bus_bin=["bus", "bin"])

            supcar = car.split("-", 2)[0]
            if supcar == "offwind":
                distance = ds["average_distance"].to_pandas()
                distance.index = distance.index.map(flatten)
                submarine_cost = costs.at[car + "-connection-submarine", "capital_cost"]
                underground_cost = costs.at[
                    car + "-connection-underground", "capital_cost"
                ]
                connection_cost = line_length_factor * (
                    distance * submarine_cost + landfall_length * underground_cost
                )

                capital_cost = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[car + "-station", "capital_cost"]
                    + connection_cost
                )
                logger.info(
                    f"Added connection cost of {connection_cost.min():0.0f}-{connection_cost.max():0.0f} Eur/MW/a to {car}"
                )
            else:
                capital_cost = costs.at[car, "capital_cost"]

            buses = ds.indexes["bus_bin"].get_level_values("bus")
            bus_bins = ds.indexes["bus_bin"].map(flatten)

            p_nom_max = ds["p_nom_max"].to_pandas()
            p_nom_max.index = p_nom_max.index.map(flatten)

            p_max_pu = ds["profile"].to_pandas()
            p_max_pu.columns = p_max_pu.columns.map(flatten)

            if not ppl.query("carrier == @car").empty:
                caps = ppl.query("carrier == @car").groupby("bus").p_nom.sum()
                caps = pd.Series(data=caps, index=ds.indexes["bus"]).fillna(0)
            else:
                caps = pd.Series(index=ds.indexes["bus"]).fillna(0)
            caps.index = caps.index.map(flatten)

            n.add(
                "Generator",
                bus_bins,
                suffix=" " + car,
                bus=buses,
                carrier=car,
                p_nom=caps,
                p_nom_min=caps,
                p_nom_extendable=car in extendable_carriers["Generator"],
                p_nom_max=p_nom_max,
                marginal_cost=costs.at[supcar, "marginal_cost"],
                capital_cost=capital_cost,
                efficiency=costs.at[supcar, "efficiency"],
                p_max_pu=p_max_pu,
                lifetime=costs.at[supcar, "lifetime"],
            )


def attach_conventional_generators(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    conventional_carriers: list,
    extendable_carriers: dict,
    conventional_params: dict,
    conventional_inputs: dict,
    unit_commitment: pd.DataFrame = None,
    fuel_price: pd.DataFrame = None,
):
    """
    Attach conventional generators to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the generators to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    ppl : pd.DataFrame
        DataFrame containing the power plant data.
    conventional_carriers : list
        List of conventional energy carriers.
    extendable_carriers : dict
        Dictionary of extendable energy carriers.
    conventional_params : dict
        Dictionary of conventional generator parameters.
    conventional_inputs : dict
        Dictionary of conventional generator inputs.
    unit_commitment : pd.DataFrame, optional
        DataFrame containing unit commitment data, by default None.
    fuel_price : pd.DataFrame, optional
        DataFrame containing fuel price data, by default None.
    """
    carriers = list(set(conventional_carriers) | set(extendable_carriers["Generator"]))

    ppl = ppl.query("carrier in @carriers")

    # reduce carriers to those in power plant dataset
    carriers = list(set(carriers) & set(ppl.carrier.unique()))
    add_missing_carriers(n, carriers)
    add_co2_emissions(n, costs, carriers)

    if unit_commitment is not None:
        committable_attrs = ppl.carrier.isin(unit_commitment).to_frame("committable")
        for attr in unit_commitment.index:
            default = n.component_attrs["Generator"].loc[attr, "default"]
            committable_attrs[attr] = ppl.carrier.map(unit_commitment.loc[attr]).fillna(
                default
            )
    else:
        committable_attrs = {}

    if fuel_price is not None:
        fuel_price = fuel_price.assign(
            OCGT=fuel_price["gas"], CCGT=fuel_price["gas"]
        ).drop("gas", axis=1)
        missing_carriers = list(set(carriers) - set(fuel_price))
        fuel_price = fuel_price.assign(**costs.fuel[missing_carriers])
        fuel_price = fuel_price.reindex(ppl.carrier, axis=1)
        fuel_price.columns = ppl.index
        marginal_cost = fuel_price.div(ppl.efficiency).add(ppl.carrier.map(costs.VOM))
    else:
        marginal_cost = ppl.marginal_cost

    # Define generators using modified ppl DataFrame
    caps = ppl.groupby("carrier").p_nom.sum().div(1e3).round(2)
    logger.info(f"Adding {len(ppl)} generators with capacities [GW]pp \n{caps}")

    n.add(
        "Generator",
        ppl.index,
        carrier=ppl.carrier,
        bus=ppl.bus,
        p_nom_min=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom_extendable=ppl.carrier.isin(extendable_carriers["Generator"]),
        efficiency=ppl.efficiency,
        marginal_cost=marginal_cost,
        capital_cost=ppl.capital_cost,
        build_year=ppl.build_year,
        lifetime=ppl.lifetime,
        **committable_attrs,
    )

    for carrier in set(conventional_params) & set(carriers):
        # Generators with technology affected
        idx = n.generators.query("carrier == @carrier").index

        for attr in list(set(conventional_params[carrier]) & set(n.generators)):
            values = conventional_params[carrier][attr]

            if f"conventional_{carrier}_{attr}" in conventional_inputs:
                # Values affecting generators of technology k country-specific
                # First map generator buses to countries; then map countries to p_max_pu
                values = pd.read_csv(
                    conventional_inputs[f"conventional_{carrier}_{attr}"], index_col=0
                ).iloc[:, 0]
                bus_values = n.buses.country.map(values)
                n.generators.update(
                    {attr: n.generators.loc[idx].bus.map(bus_values).dropna()}
                )
            else:
                # Single value affecting all generators of technology k indiscriminantely of country
                n.generators.loc[idx, attr] = values


def attach_hydro(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    profile_hydro: str,
    hydro_capacities: str,
    carriers: list,
    **params,
):
    """
    Attach hydro generators and storage units to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the hydro units to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    ppl : pd.DataFrame
        DataFrame containing the power plant data.
    profile_hydro : str
        Path to the hydro profile data.
    hydro_capacities : str
        Path to the hydro capacities data.
    carriers : list
        List of hydro energy carriers.
    **params :
        Additional parameters for hydro units.
    """
    add_missing_carriers(n, carriers)
    add_co2_emissions(n, costs, carriers)

    ror = ppl.query('carrier == "ror"')
    phs = ppl.query('carrier == "PHS"')
    hydro = ppl.query('carrier == "hydro"')

    country = ppl["bus"].map(n.buses.country).rename("country")

    inflow_idx = ror.index.union(hydro.index)
    if not inflow_idx.empty:
        dist_key = ppl.loc[inflow_idx, "p_nom"].groupby(country).transform(normed)

        with xr.open_dataarray(profile_hydro) as inflow:
            inflow_countries = pd.Index(country[inflow_idx])
            missing_c = inflow_countries.unique().difference(
                inflow.indexes["countries"]
            )
            assert missing_c.empty, (
                f"'{profile_hydro}' is missing "
                f"inflow time-series for at least one country: {', '.join(missing_c)}"
            )

            inflow_t = (
                inflow.sel(countries=inflow_countries)
                .rename({"countries": "name"})
                .assign_coords(name=inflow_idx)
                .transpose("time", "name")
                .to_pandas()
                .multiply(dist_key, axis=1)
            )

    if "ror" in carriers and not ror.empty:
        n.add(
            "Generator",
            ror.index,
            carrier="ror",
            bus=ror["bus"],
            p_nom=ror["p_nom"],
            efficiency=costs.at["ror", "efficiency"],
            capital_cost=costs.at["ror", "capital_cost"],
            weight=ror["p_nom"],
            p_max_pu=(
                inflow_t[ror.index]  # pylint: disable=E0606
                .divide(ror["p_nom"], axis=1)
                .where(lambda df: df <= 1.0, other=1.0)
            ),
        )

    if "PHS" in carriers and not phs.empty:
        # fill missing max hours to params value and
        # assume no natural inflow due to lack of data
        max_hours = params.get("PHS_max_hours", 6)
        phs = phs.replace({"max_hours": {0: max_hours, np.nan: max_hours}})
        n.add(
            "StorageUnit",
            phs.index,
            carrier="PHS",
            bus=phs["bus"],
            p_nom=phs["p_nom"],
            capital_cost=costs.at["PHS", "capital_cost"],
            max_hours=phs["max_hours"],
            efficiency_store=np.sqrt(costs.at["PHS", "efficiency"]),
            efficiency_dispatch=np.sqrt(costs.at["PHS", "efficiency"]),
            cyclic_state_of_charge=True,
        )

    if "hydro" in carriers and not hydro.empty:
        hydro_max_hours = params.get("hydro_max_hours")

        assert hydro_capacities is not None, "No path for hydro capacities given."

        hydro_stats = pd.read_csv(
            hydro_capacities, comment="#", na_values="-", index_col=0
        )
        e_target = hydro_stats["E_store[TWh]"].clip(lower=0.2) * 1e6
        e_installed = hydro.eval("p_nom * max_hours").groupby(hydro.country).sum()
        e_missing = e_target - e_installed
        missing_mh_i = hydro.query("max_hours.isnull() or max_hours == 0").index
        # some countries may have missing storage capacity but only one plant
        # which needs to be scaled to the target storage capacity
        missing_mh_single_i = hydro.index[
            ~hydro.country.duplicated() & hydro.country.isin(e_missing.dropna().index)
        ]
        missing_mh_i = missing_mh_i.union(missing_mh_single_i)

        if hydro_max_hours == "energy_capacity_totals_by_country":
            # watch out some p_nom values like IE's are totally underrepresented
            max_hours_country = (
                e_missing / hydro.loc[missing_mh_i].groupby("country").p_nom.sum()
            )

        elif hydro_max_hours == "estimate_by_large_installations":
            max_hours_country = (
                hydro_stats["E_store[TWh]"] * 1e3 / hydro_stats["p_nom_discharge[GW]"]
            )
        else:
            raise ValueError(f"Unknown hydro_max_hours method: {hydro_max_hours}")

        max_hours_country.clip(0, inplace=True)

        missing_countries = pd.Index(hydro["country"].unique()).difference(
            max_hours_country.dropna().index
        )
        if not missing_countries.empty:
            logger.warning(
                f"Assuming max_hours=6 for hydro reservoirs in the countries: {', '.join(missing_countries)}"
            )
        hydro_max_hours = hydro.max_hours.where(
            (hydro.max_hours > 0) & ~hydro.index.isin(missing_mh_single_i),
            hydro.country.map(max_hours_country),
        ).fillna(6)

        if params.get("flatten_dispatch", False):
            buffer = params.get("flatten_dispatch_buffer", 0.2)
            average_capacity_factor = inflow_t[hydro.index].mean() / hydro["p_nom"]
            p_max_pu = (average_capacity_factor + buffer).clip(upper=1)
        else:
            p_max_pu = 1

        n.add(
            "StorageUnit",
            hydro.index,
            carrier="hydro",
            bus=hydro["bus"],
            p_nom=hydro["p_nom"],
            max_hours=hydro_max_hours,
            capital_cost=costs.at["hydro", "capital_cost"],
            marginal_cost=costs.at["hydro", "marginal_cost"],
            p_max_pu=p_max_pu,  # dispatch
            p_min_pu=0.0,  # store
            efficiency_dispatch=costs.at["hydro", "efficiency"],
            efficiency_store=0.0,
            cyclic_state_of_charge=True,
            inflow=inflow_t.loc[:, hydro.index],
        )


def attach_GEM_renewables(
    n: pypsa.Network, tech_map: dict[str, list[str]], smk_inputs: list[str]
) -> None:
    """
    Attach renewable capacities from the GEM dataset to the network.

    Args:
    - n: The PyPSA network to attach the capacities to.
    - tech_map: A dictionary mapping fuel types to carrier names.

    Returns:
    - None
    """
    tech_string = ", ".join(tech_map.values())
    logger.info(f"Using GEM renewable capacities for carriers {tech_string}.")

    df = pm.data.GEM().powerplant.convert_country_to_alpha2()
    technology_b = ~df.Technology.isin(["Onshore", "Offshore"])
    df["Fueltype"] = df.Fueltype.where(technology_b, df.Technology).replace(
        {"Solar": "PV"}
    )

    for fueltype, carrier in tech_map.items():
        fn = smk_inputs.get(f"class_regions_{carrier}")
        class_regions = gpd.read_file(fn)

        df_fueltype = df.query("Fueltype == @fueltype")
        geometry = gpd.points_from_xy(df_fueltype.lon, df_fueltype.lat)
        caps = gpd.GeoDataFrame(df_fueltype, geometry=geometry, crs=4326)
        caps = caps.sjoin(class_regions)
        caps = caps.groupby(["bus", "bin"]).Capacity.sum()
        caps.index = caps.index.map(flatten) + " " + carrier

        n.generators.update({"p_nom": caps.dropna()})
        n.generators.update({"p_nom_min": caps.dropna()})


def estimate_renewable_capacities(
    n: pypsa.Network,
    year: int,
    tech_map: dict,
    expansion_limit: bool,
    countries: list,
):
    """
    Estimate a different between renewable capacities in the network and
    reported country totals from IRENASTAT dataset. Distribute the difference
    with a heuristic.

    Heuristic: n.generators_t.p_max_pu.mean() * n.generators.p_nom_max

    Args:
    - n: The PyPSA network.
    - year: The year of optimisation.
    - tech_map: A dictionary mapping fuel types to carrier names.
    - expansion_limit: Boolean value from config file
    - countries: A list of country codes to estimate capacities for.

    Returns:
    - None
    """
    if not len(countries) or not len(tech_map):
        return

    capacities = pm.data.IRENASTAT().powerplant.convert_country_to_alpha2()
    capacities = capacities.query(
        "Year == @year and Technology in @tech_map and Country in @countries"
    )
    capacities = capacities.groupby(["Technology", "Country"]).Capacity.sum()

    logger.info(
        f"Heuristics applied to distribute renewable capacities [GW]: "
        f"\n{capacities.groupby('Technology').sum().div(1e3).round(2)}"
    )

    for ppm_technology, tech in tech_map.items():
        tech_i = n.generators.query("carrier == @tech").index
        if ppm_technology in capacities.index.get_level_values("Technology"):
            stats = capacities.loc[ppm_technology].reindex(countries, fill_value=0.0)
        else:
            stats = pd.Series(0.0, index=countries)
        country = n.generators.bus[tech_i].map(n.buses.country)
        existent = n.generators.p_nom[tech_i].groupby(country).sum()
        missing = stats - existent
        dist = n.generators_t.p_max_pu.mean() * n.generators.p_nom_max

        n.generators.loc[tech_i, "p_nom"] += (
            dist[tech_i]
            .groupby(country)
            .transform(lambda s: normed(s) * missing[s.name])
            .where(lambda s: s > 0.1, 0.0)  # only capacities above 100kW
        )
        n.generators.loc[tech_i, "p_nom_min"] = n.generators.loc[tech_i, "p_nom"]

        if expansion_limit:
            assert np.isscalar(expansion_limit)
            logger.info(
                f"Reducing capacity expansion limit to {expansion_limit * 100:.2f}% of installed capacity."
            )
            n.generators.loc[tech_i, "p_nom_max"] = (
                expansion_limit * n.generators.loc[tech_i, "p_nom_min"]
            )


def attach_storageunits(
    n: pypsa.Network,
    costs: pd.DataFrame,
    extendable_carriers: dict,
    max_hours: dict,
):
    """
    Attach storage units to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the storage units to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    extendable_carriers : dict
        Dictionary of extendable energy carriers.
    max_hours : dict
        Dictionary of maximum hours for storage units.
    """
    carriers = extendable_carriers["StorageUnit"]

    n.add("Carrier", carriers)

    buses_i = n.buses.index

    lookup_store = {"H2": "electrolysis", "battery": "battery inverter"}
    lookup_dispatch = {"H2": "fuel cell", "battery": "battery inverter"}

    for carrier in carriers:
        roundtrip_correction = 0.5 if carrier == "battery" else 1

        n.add(
            "StorageUnit",
            buses_i,
            " " + carrier,
            bus=buses_i,
            carrier=carrier,
            p_nom_extendable=True,
            capital_cost=costs.at[carrier, "capital_cost"],
            marginal_cost=costs.at[carrier, "marginal_cost"],
            efficiency_store=costs.at[lookup_store[carrier], "efficiency"]
            ** roundtrip_correction,
            efficiency_dispatch=costs.at[lookup_dispatch[carrier], "efficiency"]
            ** roundtrip_correction,
            max_hours=max_hours[carrier],
            cyclic_state_of_charge=True,
        )


def attach_stores(
    n: pypsa.Network,
    costs: pd.DataFrame,
    extendable_carriers: dict,
):
    """
    Attach stores to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the stores to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    extendable_carriers : dict
        Dictionary of extendable energy carriers.
    """
    carriers = extendable_carriers["Store"]

    n.add("Carrier", carriers)

    buses_i = n.buses.index

    if "H2" in carriers:
        h2_buses_i = n.add("Bus", buses_i + " H2", carrier="H2", location=buses_i)

        n.add(
            "Store",
            h2_buses_i,
            bus=h2_buses_i,
            carrier="H2",
            e_nom_extendable=True,
            e_cyclic=True,
            capital_cost=costs.at["hydrogen storage underground", "capital_cost"],
        )

        n.add(
            "Link",
            h2_buses_i + " Electrolysis",
            bus0=buses_i,
            bus1=h2_buses_i,
            carrier="H2 electrolysis",
            p_nom_extendable=True,
            efficiency=costs.at["electrolysis", "efficiency"],
            capital_cost=costs.at["electrolysis", "capital_cost"],
            marginal_cost=costs.at["electrolysis", "marginal_cost"],
        )

        n.add(
            "Link",
            h2_buses_i + " Fuel Cell",
            bus0=h2_buses_i,
            bus1=buses_i,
            carrier="H2 fuel cell",
            p_nom_extendable=True,
            efficiency=costs.at["fuel cell", "efficiency"],
            # NB: fixed cost is per MWel
            capital_cost=costs.at["fuel cell", "capital_cost"]
            * costs.at["fuel cell", "efficiency"],
            marginal_cost=costs.at["fuel cell", "marginal_cost"],
        )

    if "battery" in carriers:
        b_buses_i = n.add(
            "Bus", buses_i + " battery", carrier="battery", location=buses_i
        )

        n.add(
            "Store",
            b_buses_i,
            bus=b_buses_i,
            carrier="battery",
            e_cyclic=True,
            e_nom_extendable=True,
            capital_cost=costs.at["battery storage", "capital_cost"],
            marginal_cost=costs.at["battery", "marginal_cost"],
        )

        n.add("Carrier", ["battery charger", "battery discharger"])

        n.add(
            "Link",
            b_buses_i + " charger",
            bus0=buses_i,
            bus1=b_buses_i,
            carrier="battery charger",
            # the efficiencies are "round trip efficiencies"
            efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
            capital_cost=costs.at["battery inverter", "capital_cost"],
            p_nom_extendable=True,
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        )

        n.add(
            "Link",
            b_buses_i + " discharger",
            bus0=b_buses_i,
            bus1=buses_i,
            carrier="battery discharger",
            efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
            p_nom_extendable=True,
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("add_electricity", clusters=100)
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    params = snakemake.params
    max_hours = params.electricity["max_hours"]
    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in params.renewable.items()
        if "landfall_length" in settings.keys()
    }

    n = pypsa.Network(snakemake.input.base_network)

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    n.set_snapshots(time)

    Nyears = n.snapshot_weightings.objective.sum() / 8760.0

    costs = load_costs(
        snakemake.input.tech_costs,
        params.costs,
        max_hours,
        Nyears,
    )

    ppl = load_and_aggregate_powerplants(
        snakemake.input.powerplants,
        costs,
        params.consider_efficiency_classes,
        params.aggregation_strategies,
        params.exclude_carriers,
    )

    attach_load(
        n,
        snakemake.input.load,
        snakemake.input.busmap,
        params.scaling_factor,
    )

    set_transmission_costs(
        n,
        costs,
        params.line_length_factor,
        params.link_length_factor,
    )

    renewable_carriers = set(params.electricity["renewable_carriers"])
    extendable_carriers = params.electricity["extendable_carriers"]
    conventional_carriers = params.electricity["conventional_carriers"]
    conventional_inputs = {
        k: v for k, v in snakemake.input.items() if k.startswith("conventional_")
    }

    if params.conventional["unit_commitment"]:
        unit_commitment = pd.read_csv(snakemake.input.unit_commitment, index_col=0)
    else:
        unit_commitment = None

    if params.conventional["dynamic_fuel_price"]:
        fuel_price = pd.read_csv(
            snakemake.input.fuel_price, index_col=0, header=0, parse_dates=True
        )
        fuel_price = fuel_price.reindex(n.snapshots).ffill()
    else:
        fuel_price = None

    attach_conventional_generators(
        n,
        costs,
        ppl,
        conventional_carriers,
        extendable_carriers,
        params.conventional,
        conventional_inputs,
        unit_commitment=unit_commitment,
        fuel_price=fuel_price,
    )

    attach_wind_and_solar(
        n,
        costs,
        ppl,
        snakemake.input,
        renewable_carriers,
        extendable_carriers,
        params.line_length_factor,
        landfall_lengths,
    )

    if "hydro" in renewable_carriers:
        p = params.renewable["hydro"]
        carriers = p.pop("carriers", [])
        attach_hydro(
            n,
            costs,
            ppl,
            snakemake.input.profile_hydro,
            snakemake.input.hydro_capacities,
            carriers,
            **p,
        )

    estimate_renewable_caps = params.electricity["estimate_renewable_capacities"]
    if estimate_renewable_caps["enable"]:
        if params.foresight != "overnight":
            logger.info(
                "Skipping renewable capacity estimation because they are added later "
                "in rule `add_existing_baseyear` with foresight mode 'myopic'."
            )
        else:
            tech_map = estimate_renewable_caps["technology_mapping"]
            expansion_limit = estimate_renewable_caps["expansion_limit"]
            year = estimate_renewable_caps["year"]

            if estimate_renewable_caps["from_gem"]:
                attach_GEM_renewables(n, tech_map, snakemake.input)

            estimate_renewable_capacities(
                n, year, tech_map, expansion_limit, params.countries
            )

    update_p_nom_max(n)

    attach_storageunits(n, costs, extendable_carriers, max_hours)
    attach_stores(n, costs, extendable_carriers)

    sanitize_carriers(n, snakemake.config)
    if "location" in n.buses:
        sanitize_locations(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
