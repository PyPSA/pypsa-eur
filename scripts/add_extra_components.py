# coding: utf-8
"""
Adds extra extendable generators, store and/or storage units to the clustered and simplified network.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        year:
        USD2013_to_EUR2013:
        dicountrate:
        emission_prices:

    electricity:
        max_hours:
        marginal_cost:
        capital_cost:
        conventional_carriers:
        co2limit:
        extendable_carriers:
            Generator:
            StorageUnit:
        estimate_renewable_capacities_from_capacity_stats:

    load:
        scaling_factor:

    renewable: (keys)
        hydro:
            carriers:
            hydro_max_hours:
            hydro_capital_cost:

    lines:
        length_factor:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`, :ref:`load_cf`, :ref:`renewable_cf`, :ref:`lines_cf`

Inputs
------

- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``data/bundle/hydro_capacities.csv``: Hydropower plant store/discharge power capacities, energy storage capacity, and average hourly inflow by country.

    .. image:: ../img/hydrocapacities.png
        :scale: 34 %

- ``data/geth2015_hydro_capacities.csv``: alternative to capacities above; NOT CURRENTLY USED!
- ``data/bundle/time_series_60min_singleindex_filtered.csv``: Hourly per-country load profiles since 2010 from the `ENTSO-E statistical database <https://www.entsoe.eu/data/power-stats/hourly_load/>`_

    .. image:: ../img/load-box.png
        :scale: 33 %

    .. image:: ../img/load-ts.png
        :scale: 33 %

- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/nuts3_shapes.geojson``: confer :ref:`shapes`
- ``resources/powerplants.csv``: confer :ref:`powerplants`
- ``resources/profile_{}.nc``: all technologies in ``config["renewables"].keys()``, confer :ref:`renewableprofiles`.
- ``networks/base.nc``: confer :ref:`base`

Outputs
-------

- ``networks/elec.nc``:

    .. image:: ../img/elec.png
            :scale: 33 %

Description
-----------

The rule :mod:`add_electricity` ties all the different data inputs from the preceding rules together into a detailed PyPSA network that is stored in ``networks/elec.nc``. It includes:

- today's transmission topology and transfer capacities (optionally including lines which are under construction according to the config settings ``lines: under_construction`` and ``links: under_construction``),
- today's thermal and hydro power generation capacities (for the technologies listed in the config setting ``electricity: conventional_carriers``), and
- today's load time-series (upsampled in a top-down approach according to population and gross domestic product)

It further adds extendable ``generators`` and ``storage_units`` with **zero** capacity for

- photovoltaic, onshore and AC- as well as DC-connected offshore wind installations with today's locational, hourly wind and solar capacity factors (but **no** current capacities),
- long-term hydrogen and short-term battery storage units (if listed in the config setting ``electricity: extendable_carriers``), and
- additional open- and combined-cycle gas turbines (if ``OCGT`` and/or ``CCGT`` is listed in the config setting ``electricity: extendable_carriers``)
"""

import logging
import pandas as pd
import pypsa
from add_electricity import load_costs, normed,
                            add_nice_carrier_names,
                            _add_missing_carriers_from_costs

idx = pd.IndexSlice
logger = logging.getLogger(__name__)


def attach_storageunits(n, costs):
    elec_opts = snakemake.config['electricity']
    carriers = elec_opts['extendable_carriers']['StorageUnit']
    max_hours = elec_opts['max_hours']

    _add_missing_carriers_from_costs(n, costs, carriers)

    buses_i = n.buses.index

    for carrier in carriers:
        n.madd("StorageUnit", buses_i, ' ' + carrier,
               bus=buses_i,
               carrier=carrier,
               p_nom_extendable=True,
               capital_cost=costs.at[carrier, 'capital_cost'],
               marginal_cost=costs.at[carrier, 'marginal_cost'],
               efficiency_store=costs.at[carrier, 'efficiency'],
               efficiency_dispatch=costs.at[carrier, 'efficiency'],
               max_hours=max_hours[carrier],
               cyclic_state_of_charge=True)

def attach_stores(n, costs):
    elec_opts = snakemake.config['electricity']
    carriers = elec_opts['extendable_carriers']['Store']

    _add_missing_carriers_from_costs(n, costs, carriers)

    buses_i = n.buses.index
    bus_sub_dict = {k: n.buses[k].values for k in ['x', 'y', 'country']}

    if 'H2' in carriers:
        h2_buses_i = n.madd("Bus", buses_i + " H2", carrier="H2", **bus_sub_dict)

        n.madd("Store", h2_buses_i,
               bus=h2_buses_i,
               e_nom_extendable=True,
               e_cyclic=True,
               capital_cost=costs.at["hydrogen storage", "capital_cost"])

        n.madd("Link", h2_buses_i + " Electrolysis",
               bus0=buses_i,
               bus1=h2_buses_i,
               p_nom_extendable=True,
               efficiency=costs.at["electrolysis", "efficiency"],
               capital_cost=costs.at["electrolysis", "capital_cost"])

        n.madd("Link", h2_buses_i + " Fuel Cell",
               bus0=h2_buses_i,
               bus1=buses_i,
               p_nom_extendable=True,
               efficiency=costs.at["fuel cell", "efficiency"],
               #NB: fixed cost is per MWel
               capital_cost=costs.at["fuel cell", "capital_cost"] * costs.at["fuel cell", "efficiency"])

    if 'battery' in carriers:
        b_buses_i = n.madd("Bus", buses_i + " battery", carrier="battery", **bus_sub_dict)

        n.madd("Store", b_buses_i,
               bus=b_buses_i,
               e_cyclic=True,
               e_nom_extendable=True,
               capital_cost=costs.at['battery storage', 'capital_cost'])

        n.madd("Link", b_buses_i + " charger",
               bus0=buses_i,
               bus1=b_buses_i,
               efficiency=costs.at['battery inverter', 'efficiency'],
               capital_cost=costs.at['battery inverter', 'capital_cost'],
               p_nom_extendable=True)

        n.madd("Link", b_buses_i + " discharger",
               bus0=b_buses_i,
               bus1=buses_i,
               efficiency=costs.at['battery inverter','efficiency'],
               capital_cost=costs.at['battery inverter', 'capital_cost'],
               p_nom_extendable=True)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake, Dict

        snakemake = MockSnakemake(output=['networks/elec_s_5_ec.nc'])
        snakemake.input = snakemake.expand(
            Dict(network='networks/elec_s_5.nc',
                 tech_costs='data/costs.csv'))

    logging.basicConfig(level=snakemake.config['logging_level'])

    n = pypsa.Network(snakemake.input.network)
    Nyears = n.snapshot_weightings.sum()/8760.
    costs = load_costs(Nyears, tech_costs=snakemake.input.tech_costs,
                       config=snakemake.config['costs'],
                       elec_config=snakemake.config['electricity'])

    attach_storageunits(n, costs)
    attach_stores(n, costs)

    add_nice_carrier_names(n, config=snakemake.config)

    n.export_to_netcdf(snakemake.output[0])
