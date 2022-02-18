.. _myopic:

##########################################
Myopic transition path
##########################################

The myopic code can be used to investigate progressive changes in a network, for instance, those taking place throughout a transition path. The capacities installed in a certain time step are maintained in the network until their operational lifetime expires.

The myopic approach was initially developed and used in the paper `Early decarbonisation of the European Energy system pays off (2020) <https://www.nature.com/articles/s41467-020-20015-4>`__ but the current implementation complies with the pypsa-eur-sec standard working flow and is compatible with using the higher resolution electricity transmission model `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`__ rather than a one-node-per-country model.

The current code applies the myopic approach to generators, storage technologies and links in the power sector and the space and water heating sector.

The transport sector and industry are not affected by the myopic code. In essence, the electrification of road and rail transport, the percentage of electric vehicles that allow demand-side management and vehicle-to-grid services, and the transformation in the different industrial subsectors do not evolve with time. They are kept fixed at the values specified in the configuration file. Including the transport sector and industry in the myopic code is planned for the near future.

See also other `outstanding issues <https://github.com/PyPSA/pypsa-eur-sec/issues/19#issuecomment-678194802>`_.

Configuration
=================

PyPSA-Eur-Sec has several configuration options which are collected in a config.yaml file located in the root directory. For myopic optimization, users should copy the provided default configuration ``config.default.yaml`` and make their own modifications and assumptions in the user-specific configuration file (``config.yaml``).

The following options included in the config.yaml file  are relevant for the myopic code.

To activate the myopic option select ``foresight: 'myopic'`` in ``config.yaml``.

To set the investment years which are sequentially simulated for the myopic investment planning, select for example ``planning_horizons : [2020, 2030, 2040, 2050]`` in ``config.yaml``.



**existing capacities**

Grouping years indicates the bins limits for grouping the existing capacities of different technologies

grouping_years: [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019]




**threshold capacity**

if for a technology, node, and grouping bin, the capacity is lower than threshold_capacity, it is ignored

threshold_capacity: 10




**conventional carriers**

conventional carriers indicate carriers used in the existing conventional technologies

conventional_carriers: ['lignite', 'coal', 'oil', 'uranium']



Wildcards
==============================

**{planning_horizons} wildcard**

The {planning_horizons} wildcard indicates the timesteps in which the network is optimized, e.g. planning_horizons: [2020, 2030, 2040, 2050]


Options
=============
The total carbon budget for the entire transition path can be indicated in the ``scenario.sector_opts`` in ``config.yaml``.
The carbon budget can be split among the ``planning_horizons`` following an exponential or beta decay. 
E.g. ``'cb40ex0'`` splits the a carbon budget equal to 40 GtCO_2 following an exponential decay whose initial linear growth rate $r$ is zero

$e(t) = e_0 (1+ (r+m)t) e^(-mt)$

See details in Supplementary Note 1 of the paper `Early decarbonisation of the European Energy system pays off (2020) <https://www.nature.com/articles/s41467-020-20015-4>`__

Rules overview
=================

General myopic code structure
===============================

The myopic code solves the network for the time steps included in ``planning_horizons`` in a recursive loop, so that:

1.The existing capacities (those installed before the base year are added as fixed capacities with p_nom=value, p_nom_extendable=False). E.g. for baseyear=2020, capacities installed before 2020 are added. In addition, the network comprises additional generator, storage, and link capacities with p_nom_extendable=True. The non-solved network is saved in ``results/run_name/networks/prenetworks-brownfield``.

The base year is the first element in ``planning_horizons``. Step 1 is implemented with the rule add_baseyear for the base year and with the rule add_brownfield for the remaining planning_horizons.

2.The 2020 network is optimized. The solved network is saved in ``results/run_name/networks/postnetworks``

3.For the next planning horizon, e.g. 2030, the capacities from a previous time step are added if they are still in operation (i.e., if they fulfil planning horizon <= commissioned year + lifetime). In addition, the network comprises additional generator, storage, and link capacities with p_nom_extendable=True. The non-solved network is saved in ``results/run_name/networks/prenetworks-brownfield``.

Steps 2 and 3 are solved recursively for all the planning_horizons included in ``config.yaml``.


rule add_existing baseyear
==========================

The rule add_existing_baseyear loads the network in ‘results/run_name/networks/prenetworks’ and performs the following operations:

1.Add the conventional, wind and solar power generators that were installed before the base year.

2.Add the heating capacities that were installed before the base year.

The existing conventional generators are retrieved from the `powerplants.csv file <https://pypsa-eur.readthedocs.io/en/latest/preparation/build_powerplants.html?highlight=powerplants>`__ generated by pypsa-eur which, in turn, is based on the `powerplantmatching <https://github.com/FRESNA/powerplantmatching>`__ database.

Existing wind and solar capacities are retrieved from `IRENA annual statistics <https://www.irena.org/Statistics/Download-Data>`__ and distributed among the nodes in a country proportional to capacity factor. (This will be updated to include capacity distributions closer to reality.)

Existing heating capacities are retrieved from the report `Mapping and analyses of the current and future (2020 - 2030) heating/cooling fuel deployment (fossil/renewables)
<https://ec.europa.eu/energy/studies/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment_en?redir=1>`__

The heating capacities are assumed to have a lifetime indicated by the parameter lifetime in the configuration file, e.g 25 years. They are assumed to be decommissioned linearly starting on the base year, e.g., from 2020 to 2045.

Then, the resulting network is saved in ``results/run_name/networks/prenetworks-brownfield``.

rule add_brownfield
===================

The rule add_brownfield loads the network in ``results/run_name/networks/prenetworks`` and performs the following operation:

1.Read the capacities optimized in the previous time step and add them to the network if they are still in operation (i.e., if they fulfill planning horizon < commissioned year + lifetime)

Then, the resulting network is saved in ``results/run_name/networks/prenetworks_brownfield``.
