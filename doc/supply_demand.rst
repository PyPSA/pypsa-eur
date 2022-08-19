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
* ``residential/services rural``: heating for individual buildings in rural areas, agriculture heat uses


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
=============================

Hydrogen is consumed in the industry sector (link to industry) to produce ammonia [link to ammonia industry section] and direct reduced iron (DRI) [link to DRI industry section]. Hydrogen is also consumed to produce synthetic methane [link to section “Methane supply”] and liquid hydrocarbons [link to fossil-oil based supply”] which have multiple uses in industry and other sectors. 
Hydrogen is also used for transport applications (link to transport), where it is exogenously fixed. It is used in `heavy-duty land transport <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L181>`_ and as liquified hydrogen in the shipping sector [add link to shipping sector]. Furthermore, stationary fuel cells may re-electrify hydrogen (with waste heat as a byproduct) to balance renewable fluctuations [Add a link to the section where we describe the Electricity sector and how storage is modelled there]. The waste heat from the stationary fuel cells can be used in `district-heating systems <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L256>`_.

Hydrogen supply
=============================

Today, most of the H2 consumed globally is produced from natural gas by steam methane reforming (SMR)

$$
CH_4 + H_2O → CO + 3H_2
$$

combined with a water-gas shift reaction

$$
CO + H_2O → CO_2 + H_2
$$

SMR is included `here <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L245>`_. 
PyPSA-Eur-Sec allows this route of H2 production with and without [carbon capture (CC)] (Link to section on Carbon Capture Storage and Utilization). These routes are often referred to as blue and grey hydrogen. Here, methane input can be both of fossil or synthetic origin.

Green hydrogen can be produced by electrolysis to split water into hydrogen and oxygen

$$
2H_2O → 2H_2 + O_2
$$

For the electrolysis, alkaline electrolysers are chosen since they have lower cost and higher cumulative installed capacity than polymer electrolyte membrane (PEM) electrolysers. The techno-economic assumptions are taken from the technology-data repository. Waste heat from electrolysis is not leveraged in the model.

*Transport*

Hydrogen is transported by pipelines. H2 pipelines are endogenously generated, either via a  greenfield H2 network, or by `retrofitting natural gas pipelines <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L262>`_). Retrofitting is implemented in such a way that for every unit of decommissioned gas pipeline, a share (60% is used in [link to H2 backbone study]) of its nominal capacity (exogenously determined in the `config file <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L266>`_.) is available for hydrogen transport. When the gas network is not resolved, this input denotes the potential for gas pipelines repurposed into hydrogen pipelines.
New pipelines can be built additionally on all routes where there currently is a gas or electricity network connection. These new pipelines will be built where no sufficient retrofitting options are available. The capacities of new and repurposed pipelines are a result of the optimisation.

*Storage*

Hydrogen can be stored in overground steel tanks or `underground salt caverns <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L250>`_. For the latter, energy storage capacities in every country are limited to the potential estimation for onshore salt caverns within `50 km <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L251>`_ of shore to avoid environmental issues associated with brine solution disposal. Underground storage potentials for hydrogen in European salt caverns is acquired from `Caglayan et al. <https://doi.org/10.1016/j.ijhydene.2019.12.161>`_

Methane demand  
====================================
                  
Methane is used in individual and large-scale gas boilers, in CHP plants with and without carbon capture, in OCGT and CCGT power plants, and in some industry subsectors for the provision of high temperature heat[LINK TO INDUSTRY OVERVIEW]. Methane is not used in the trans- port sector because of engine slippage.

Methane supply                          
===================================

In addition to methane from fossil origins, the model also considers biogenic and synthetic sources. `The gas network can either be modeled, or it can be assumed that gas transport is not limited <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L261>`_. If gas infrastructure is regionally resolved, fossil gas can enter the system only at existing and planned LNG terminals, pipeline entry-points, and intra- European gas extraction sites, which are retrieved from the SciGRID Gas IGGIELGN dataset and the GEM Wiki. 
Biogas can be upgraded to methane. 
Synthetic methane can be produced by processing hydrogen and captures CO2 in the Sabatier reaction 

$$                   
CO_2 + 4H_2 → CH_4 + 2H_2O
$$

Direct power-to-methane conversion with efficient heat integration developed in the HELMETH project is also an option. The share of synthetic, biogenic and fossil methane is an optimisation result depending on the techno-economic assumptions.
                
*Methane transport*
            
The existing European gas transmission network is represented based on the SciGRID Gas IGGIELGN dataset. This dataset is based on compiled and merged data from the ENTSOG maps  and other publicly available data sources. It includes data on the capacity, diameter, pressure, length, and directionality of pipelines. Missing capacity data is conservatively inferred from the pipe diameter following conversion factors derived from an EHB report. The gas network is clustered to the selected number of model regions. Gas pipelines can be endogenously expanded or repurposed for hydrogen transport. Gas flows are represented by a lossless transport model. Methane is assumed to be transmitted without cost or capacity constraints because future demand is predicted to be low compared to available transport capacities.

The following figure shows the unclustered European gas transmission network based on the SciGRID Gas IGGIELGN dataset. Pipelines are color-coded by estimated capacities. Markers indicate entry-points, sites of fossil resource extraction, and LNG terminals.            
    
.. image:: ../graphics/gas_pipeline_figure.png


Biomass
============

Biomass supply
---------------
Biomass supply potentials for each European country are taken from the `JRC ENSPRESO database <http://data.europa.eu/89h/74ed5a04-7d74-4807-9eab-b94774309d9f>`_ where data is available for various years (2010, 2020, 2030, 2040 and 2050) and scenarios (low, medium, high). No biomass import from outside Europe is assumed. More information on the data set can be found `here <https://publications.jrc.ec.europa.eu/repository/handle/JRC98626>`_.

Solid biomass demand
=====================


The desired scenario can be selected in the pypsa-eur-sec `configuration <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L108>`_. The script for building the biomass potentials from the JREC ENSPRESO data base is located `here <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/scripts/build_biomass_potentials.py#L43>`_. Consult the script to see the keywords that specify the scenario options.

The `configuration <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L108>`_ also allows the user to define how the various types of biomass are used in the model by using the categories : biogas, solid biomass, and not included.
Feedstocks categorized as biogas, typically manure and sludge waste, are available to the model as biogas (that is upgraded to biomethane). More details below.

Feedstocks categorized as solid biomass, e.g. secondary forest residues or municipal waste can be used directly or converted to gas or liquid fuels. More details below.

Feedstocks labeled as not included are ignored by the model.
A `typical use case for biomass <https://arxiv.org/abs/2109.09563>`_ would be the medium availability scenario for 2030 where only residues from agriculture and forestry as well as biodegradable municipal waste are considered as energy feedstocks. Fuel crops are avoided because they compete with scarce land for food production, while primary wood, as well as wood chips and pellets, are avoided because of concerns about sustainability . See the supporting materials of the `paper <https://www.sciencedirect.com/science/article/pii/S1364032117302034>`_ for more details.

Solid biomass conversion and use
----------------------------------
Solid biomass can be used directly to provide process heat up to 500 C in the industry. It can also be burnt in CHP plants and boilers associated with heating systems. These technologies are described elsewhere [link to heat and industry sections].

Solid biomass can be converted to syngas if the option is enabled in the `config file <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L274>`_. In this case the model will enable the technology BioSNG both with and without the option for carbon capture [link to technology data].
Liquefaction of solid biomass `can be enabled <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L273>`_ allowing the model to convert it into liquid hydrocarbons that can replace conventional oil products. This technology also comes with and without carbon capture [link to technology data].

Transport of solid biomass
---------------------------
The transport of solid biomass can either be assumed unlimited between countries or it can be associated with a country specific cost per MWh/km. In the config file these options are toggled `here <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L270>`_. If the option is off, use of solid biomass is transport. If it is turned on, a biomass transport network will be `created <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/scripts/prepare_sector_network.py#L1803>`_ between all nodes. This network resembles road transport of biomass and the cost of transportation is a variable cost which is proportional to distance and a country specific cost per MWh/km. The latter is `estimated <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/scripts/build_biomass_transport_costs.py>`_ from the country specific costs per ton/km used in the publication `“The JRC-EU-TIMES model. Bioenergy potentials for EU and neighbouring countries” <https://publications.jrc.ec.europa.eu/repository/handle/JRC98626>`_.

Biogas transport and use
------------------------
Biogas will be aggregated into a common European resources if a gas network is not modeled explicitly, i.e., the `gas_network <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L261>`_ option is set to false. If, on the other hand, a gas network is included, the biogas potential will be associated with each node of origin.
The model can only use biogas by first upgrading it to natural gas quality [link to tech description] (bio methane) which is fed into the general gas network.



Oil-based products demand
========================
Naphtha is used as a feedstock in the chemicals industry[LINK TO CHEMICAL INDUSTRY]. Furthermore, kerosene is used as transport fuel in the aviation sector[LINK TO AVIATION SECTOR]. Non-electrified agriculture machinery also consumes gasoline. 
Land transport [LINK TO LAND TRANSPORT]  that is not electrified or converted into using H2-fuel cells also consumes oil-based products. While there is regional distribution of demand, the carrier is copperplated in the model, which means that transport costs and constraints are neglected. 


Oil-based products supply 
======================== 							
Oil-based products can be either of fossil origin or synthetically produced by combining H2 [link to hydrogen] and captured CO2 [link to carbon capture] in Fischer-Tropsch plants		

$$
𝑛CO+(2𝑛+1)H_2 → C_{n}H_{2n + 2}  +𝑛H_2O
$$

with costs as included from the `technology-data repository <https://github.com/PyPSA/technology-data/blob/master/latex_tables/tables_in_latex.pdf>`_. The waste heat from the Fischer-Tropsch process is supplied to `district heating networks <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L255>`_. The share of fossil and synthetic oil is an optimisation result depending on the techno-economic assumptions.


Oil-based transport
======================== 		
Liquid hydrocarbons are assumed to be transported freely among the model region since future demand is predicted to be low, transport costs for liquids are low and no bottlenecks are expected. 


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

Transportation
=========================
Annual energy demands for land transport, aviation and shipping for every country are retrieved from `JRC-IDEES data set <http://data.europa.eu/89h/jrc-10110-10001>`_. Below, the details of how each of these categories are treated is explained.

Land transport
-----------------

Aviation
-----------------
The `demand for aviation <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/scripts/prepare_sector_network.py#L2193>`_ includes international and domestic use. It is modeled as an oil demand since aviation consumes kerosene. This can be produced synthetically or have fossil-origin [link to oil product].

Shipping
----------------
Shipping energy demand is covered by a combination of oil and hydrogen. Other fuel options, like methanol or ammonia, are currently not included in PyPSA-Eur-Sec.The share of shipping that is assumed to be supplied by hydrogen can be selected in the config file <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L198>`_.

To estimate the `hydrogen demand <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/scripts/prepare_sector_network.py#L2090>`_, the average fuel efficiency of the fleet is used in combination with the efficiency of the fuel cell defined in the technology-data repository. The average fuel efficiency is set in the `config file <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L196>`_.

The consumed hydrogen comes from the general hydrogen bus where it can be produced by SMR, SMR+CC or electrolysers [link to hydrogen]. The fraction that is not converted into hydrogen use oil products, i.e., is connected to the general oil bus.

The energy demand for liquefaction of the hydrogen used for shipping can be `included  <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/config.default.yaml#L197>`_. If this option is selected, liquifaction will happen at the `node where the shipping demand occurs <https://github.com/PyPSA/pypsa-eur-sec/blob/3daff49c9999ba7ca7534df4e587e1d516044fc3/scripts/prepare_sector_network.py#L2064>`_.



Carbon dioxide capture, usage and sequestration (CCU/S)
=========================================================

Carbon dioxide can be captured from industry process emissions,
emissions related to industry process heat, combined heat and power
plants, and directly from the air (DAC).

Carbon dioxide can be used as an input for methanation and
Fischer-Tropsch fuels, or it can be sequestered underground.
