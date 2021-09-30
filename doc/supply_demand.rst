.. _supply_demand:

##########################################
Supply and demand
##########################################

An initial orientation to the supply and demand options in the model
PyPSA-Eur-Sec can be found in the description of the model
PyPSA-Eur-Sec-30 in the paper `Synergies of sector coupling and
transmission reinforcement in a cost-optimised, highly renewable
European energy system <https://arxiv.org/abs/1801.05290>`_ (2018).
The latest version of PyPSA-Eur-Sec differs by including biomass,
industry, industrial feedstocks, aviation, shipping, better carbon
management, carbon capture and usage/sequestration, and gas networks.

The basic supply (left column) and demand (right column) options in the model are described in this figure:

.. image:: ../graphics/multisector_figure.png



Electricity supply and demand
=============================

Electricity supply and demand follows the electricity generation and
transmission model `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_,
except that hydrogen storage is integrated into the hydrogen supply,
demand and network, and PyPSA-Eur-Sec includes CHPs.

Unlike PyPSA-Eur, PyPSA-Eur-Sec does not distribution electricity demand for industry according to population and GDP, but uses the
geographical data from the `Hotmaps Industrial Database
<https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`_.

Also unlike PyPSA-Eur, PyPSA-Eur-Sec subtracts existing electrified heating from the existing electricity demand, so that power-to-heat can be optimised separately.

The remaining electricity demand for households and services is distributed inside each country proportional to GDP and population.


Heat demand
=============================

Heat demand is split into:

* ``urban central``: large-scale district heating networks in urban areas with dense heat demand
* ``residential/services urban decentral``: heating for individual buildings in urban areas
* ``residential/services rural``: heating for individual buildings in rural areas


Heat supply
=======================

Oil and gas boilers
--------------------

Heat pumps
-------------

Either air-to-water or ground-to-water heat pumps are implemented.

They have coefficient of performance (COP) based on either the
external air or the soil hourly temperature.

Ground-source heat pumps are only allowed in rural areas because of
space constraints.

Only air-source heat pumps are allowed in urban areas. This is a
conservative assumption, since there are many possible sources of
low-temperature heat that could be tapped in cities (waste water,
rivers, lakes, seas, etc.).

Resistive heaters
--------------------


Large Combined Heat and Power (CHP) plants
--------------------------------------------

A good summary of CHP options that can be implemented in PyPSA can be found in the paper `Cost sensitivity of optimal sector-coupled district heating production systems <https://doi.org/10.1016/j.energy.2018.10.044>`_.

PyPSA-Eur-Sec includes CHP plants fuelled by methane, hydrogen and solid biomass from waste and residues.

Hydrogen CHPs are fuel cells.

Methane and biomass CHPs are based on back pressure plants operating with a fixed ratio of electricity to heat output. The methane CHP is modelled on the Danish Energy Agency (DEA) "Gas turbine simple cycle (large)" while the solid biomass CHP is based on the DEA's "09b Wood Pellets Medium".

The efficiencies of each are given on the back pressure line, where the back pressure coefficient ``c_b`` is the electricity output divided by the heat output. The plants are not allowed to deviate from the back pressure line and are implement as ``Link`` objects with a fixed ratio of heat to electricity output.


NB: The old PyPSA-Eur-Sec-30 model assumed an extraction plant (like the DEA coal CHP) for gas which has flexible production of heat and electricity within the feasibility diagram of Figure 4 in the `Synergies paper <https://arxiv.org/abs/1801.05290>`_. We have switched to the DEA back pressure plants since these are more common for smaller plants for biomass, and because the extraction plants were on the back pressure line for 99.5% of the time anyway. The plants were all changed to back pressure in PyPSA-Eur-Sec v0.4.0.


Micro-CHP for individual buildings
-----------------------------------

Optional.

Waste heat from Fuel Cells, Methanation and Fischer-Tropsch plants
-------------------------------------------------------------------


Solar thermal collectors
-------------------------

Thermal energy storage using hot water tanks
---------------------------------------------

Small for decentral applications.

Big water pit storage for district heating.

.. _retro:

Retrofitting of the thermal envelope of buildings
===================================================
Co-optimising building renovation is only enabled if in the ``config.yaml`` the
option :mod:`retro_endogen: True`. To reduce the computational burden
default setting is

.. literalinclude:: ../config.default.yaml
    :language: yaml
    :lines: 134-135

Renovation of the thermal envelope reduces the space heating demand and is
optimised at each node for every heat bus. Renovation measures through additional
insulation material and replacement of energy inefficient windows are considered.

In a first step, costs per energy savings are estimated in :mod:`build_retro_cost.py`.
They depend on the insulation condition of the building stock and costs for
renovation of the building elements.
In a second step, for those cost per energy savings two possible renovation
strengths are determined: a moderate renovation with lower costs and lower
maximum possible space heat savings, and an ambitious renovation with associated
higher costs and higher efficiency gains. They are added by step-wise
linearisation in form of two additional generations in
:mod:`prepare_sector_network.py`.

Settings in the config.yaml concerning the endogenously optimisation of building
renovation

.. literalinclude:: ../config.default.yaml
    :language: yaml
    :lines: 136-140

Further information are given in the publication

`Mitigating heat demand peaks in buildings in a highly renewable European energy system, (2021) <https://arxiv.org/abs/2012.01831>`_.


Hydrogen demand
==================

Stationary fuel cell CHP.

Transport applications.

Industry (ammonia, precursor to hydrocarbons for chemicals and iron/steel).


Hydrogen supply
=================

Steam Methane Reforming (SMR), SMR+CCS, electrolysers.


Methane demand
==================

Can be used in boilers, in CHPs, in industry for high temperature heat, in OCGT.

Not used in transport because of engine slippage.

Methane supply
=================

Fossil, biogas, Sabatier (hydrogen to methane), HELMETH (directly power to methane with efficient heat integration).


Solid biomass demand
=====================

Solid biomass provides process heat up to 500 Celsius in industry, as well as feeding CHP plants in district heating networks.

Solid biomass supply
=====================

Only wastes and residues from the JRC ENSPRESO biomass dataset.


Oil product demand
=====================

Transport fuels and naphtha as a feedstock for the chemicals industry.

Oil product supply
======================

Fossil or Fischer-Tropsch.


Industry demand
================

Based on materials demand from JRC-IDEES and other sources such as the USGS for ammonia.

Industry is split into many sectors, including iron and steel, ammonia, other basic chemicals, cement, non-metalic minerals, alumuninium, other non-ferrous metals, pulp, paper and printing, food, beverages and tobacco, and other more minor sectors.

Inside each country the industrial demand is distributed using the `Hotmaps Industrial Database <https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`_.


Industry supply
================

Process switching (e.g. from blast furnaces to direct reduction and electric arc furnaces for steel) is defined exogenously.

Fuel switching for process heat is mostly also done exogenously.

Solid biomass is used for up to 500 Celsius, mostly in paper and pulp and food and beverages.

Higher temperatures are met with methane.


Carbon dioxide capture, usage and sequestration (CCU/S)
=========================================================

Carbon dioxide can be captured from industry process emissions,
emissions related to industry process heat, combined heat and power
plants, and directly from the air (DAC).

Carbon dioxide can be used as an input for methanation and
Fischer-Tropsch fuels, or it can be sequestered underground.
