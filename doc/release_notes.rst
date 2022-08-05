##########################################
Release Notes
##########################################

Future release
==============

.. note::
  This unreleased version currently may require the master branches of PyPSA, PyPSA-Eur, and the technology-data repository.

This release includes the addition of the European gas transmission network and
incorporates retrofitting options to hydrogen.

**Gas Transmission Network**

* New rule ``retrieve_gas_infrastructure_data`` that downloads and extracts the
  SciGRID_gas `IGGIELGN <https://zenodo.org/record/4767098>`_ dataset from zenodo.
  It includes data on the transmission routes, pipe diameters, capacities, pressure,
  and whether the pipeline is bidirectional and carries H-Gas or L-Gas.

* New rule ``build_gas_network`` processes and cleans the pipeline data from SciGRID_gas.
  Missing or uncertain pipeline capacities can be inferred by diameter.

* New rule ``build_gas_input_locations`` compiles the LNG import capacities
  (including planned projects from gem.wiki), pipeline entry capacities and
  local production capacities for each region of the model. These are the
  regions where fossil gas can eventually enter the model.

* New rule ``cluster_gas_network`` that clusters the gas transmission network
  data to the model resolution. Cross-regional pipeline capacities are aggregated
  (while pressure and diameter compability is ignored), intra-regional pipelines
  are dropped. Lengths are recalculated based on the regions' centroids.

* With the option ``sector: gas_network:``, the existing gas network is
  added with a lossless transport model. A length-weighted `k-edge augmentation
  algorithm
  <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation>`_
  can be run to add new candidate gas pipelines such that all regions of the
  model can be connected to the gas network. The number of candidates can be
  controlled via the setting ``sector: gas_network_connectivity_upgrade:``. When
  the gas network is activated, all the gas demands are regionally disaggregated
  as well.

* New constraint allows endogenous retrofitting of gas pipelines to hydrogen pipelines.
  This option is activated via the setting ``sector: H2_retrofit:``. For every
  unit of gas pipeline capacity dismantled, ``sector:
  H2_retrofit_capacity_per_CH4`` units are made available as hydrogen pipeline
  capacity in the corresponding corridor. These repurposed hydrogen pipelines
  have lower costs than new hydrogen pipelines. Both new and repurposed pipelines
  can be built simultaneously. The retrofitting option ``sector: H2_retrofit:`` also works
  with a copperplated methane infrastructure, i.e. when ``sector: gas_network: false``.

* New hydrogen pipelines can now be built where there are already power or gas
  transmission routes. Previously, only the electricity transmission routes were
  considered.

**New features and functionality**


* Add option for biomass boilers (wood pellets) for decentral heating

* Add option for BioSNG (methane from biomass) with and without CC 

* Add option for BtL (Biomass to liquid fuel/oil) with and without CC 

* Units are assigned to the buses. These only provide a better understanding. The specifications of the units are not taken into account in the optimisation, which means that no automatic conversion of units takes place.

* Option ``retrieve_sector_databundle`` to automatically retrieve and extract data bundle.

* Add regionalised hydrogen salt cavern storage potentials from `Technical Potential of Salt Caverns for Hydrogen Storage in Europe <https://doi.org/10.20944/preprints201910.0187.v1>`_.

* Add option to sweep the global CO2 sequestration potentials with keyword ``seq200`` in the ``{sector_opts}`` wildcard (for limit of 200 Mt CO2).

* Updated `data bundle <https://zenodo.org/record/5824485/files/pypsa-eur-sec-data-bundle.tar.gz>`_ that includes the hydrogan salt cavern storage potentials.

**Bugfixes**

* The CO2 sequestration limit implemented as GlobalConstraint (introduced in the previous version)
  caused a failure to read in the shadow prices of other global constraints.


PyPSA-Eur-Sec 0.6.0 (4 October 2021)
====================================

This release includes
improvements regarding the basic chemical production,
the addition of plastics recycling,
the addition of the agriculture, forestry and fishing sector,
more regionally resolved biomass potentials,
CO2 pipeline transport and storage, and
more options in setting exogenous transition paths,
besides many performance improvements.

This release is known to work with `PyPSA-Eur
<https://github.com/PyPSA/pypsa-eur>`_ Version 0.4.0, `Technology Data
<https://github.com/PyPSA/technology-data>`_ Version 0.3.0 and
`PyPSA <https://github.com/PyPSA/PyPSA>`_ Version 0.18.0.

Please note that the data bundle has also been updated.


**General**

* With this release, we change the license from copyleft GPLv3 to the more
  liberal MIT license with the consent of all contributors.


**New features and functionality**

* Distinguish costs for home battery storage and inverter from utility-scale
  battery costs.

* Separate basic chemicals into HVC (high-value chemicals), chlorine, methanol and ammonia
  [`#166 <https://github.com/PyPSA/PyPSA-Eur-Sec/pull/166>`_].

* Add option to specify reuse, primary production, and mechanical and chemical
  recycling fraction of platics
  [`#166 <https://github.com/PyPSA/PyPSA-Eur-Sec/pull/166>`_].

* Include energy demands and CO2 emissions for the agriculture, forestry and fishing sector.
  It is included by default through the option ``A`` in the ``sector_opts`` wildcard.
  Part of the emissions (1.A.4.c) was previously assigned to "industry non-elec" in the ``co2_totals.csv``.
  Hence, excluding the agriculture sector will now lead to a tighter CO2 limit.
  Energy demands are taken from the JRC IDEES database (missing countries filled with eurostat data)
  and are split into
  electricity (lighting, ventilation, specific electricity uses, pumping devices (electric)),
  heat (specific heat uses, low enthalpy heat)
  machinery oil (motor drives, farming machine drives, pumping devices (diesel)).
  Heat demand is assigned at "services rural heat" buses.
  Electricity demands are added to low-voltage buses.
  Time series for demands are constant and distributed inside countries by population
  [`#147 <https://github.com/PyPSA/PyPSA-Eur-Sec/pull/147>`_].

* Include today's district heating shares in myopic optimisation and add option
  to specify exogenous path for district heating share increase under ``sector:
  district_heating:`` [`#149 <https://github.com/PyPSA/PyPSA-Eur-Sec/pull/149>`_].

* Added option for hydrogen liquefaction costs for hydrogen demand in shipping.
  This introduces a new ``H2 liquid`` bus at each location. It is activated via
  ``sector: shipping_hydrogen_liquefaction: true``.

* The share of shipping transformed into hydrogen fuel cell can be now defined
  for different years in the ``config.yaml`` file. The carbon emission from the
  remaining share is treated as a negative load on the atmospheric carbon dioxide
  bus, just like aviation and land transport emissions.

* The transformation of the Steel and Aluminium production can be now defined
  for different years in the ``config.yaml`` file.

* Include the option to alter the maximum energy capacity of a store via the
  ``carrier+factor`` in the ``{sector_opts}`` wildcard. This can be useful for
  sensitivity analyses. Example: ``co2 stored+e2`` multiplies the ``e_nom_max`` by
  factor 2. In this example, ``e_nom_max`` represents the CO2 sequestration
  potential in Europe.

* Use `JRC ENSPRESO database <https://data.jrc.ec.europa.eu/dataset/74ed5a04-7d74-4807-9eab-b94774309d9f>`_ to
  spatially disaggregate biomass potentials to PyPSA-Eur regions based on
  overlaps with NUTS2 regions from ENSPRESO (proportional to area) (`#151
  <https://github.com/PyPSA/pypsa-eur-sec/pull/151>`_).

* Add option to regionally disaggregate biomass potential to individual nodes
  (previously given per country, then distributed by population density within)
  and allow the transport of solid biomass. The transport costs are determined
  based on the `JRC-EU-Times Bioenergy report
  <http://dx.doi.org/10.2790/01017>`_ in the new optional rule
  ``build_biomass_transport_costs``. Biomass transport can be activated with the
  setting ``sector: biomass_transport: true``.

* Add option to regionally resolve CO2 storage and add CO2 pipeline transport
  because geological storage potential,
  CO2 utilisation sites and CO2 capture sites may be separated. The CO2 network
  is built from zero based on the topology of the electricity grid (greenfield).
  Pipelines are assumed to be bidirectional and lossless. Furthermore, neither
  retrofitting of natural gas pipelines (required pressures are too high, 80-160
  bar vs <80 bar) nor other modes of CO2 transport (by ship, road or rail) are
  considered. The regional representation of CO2 is activated with the config
  setting ``sector: co2_network: true`` but is deactivated by default. The
  global limit for CO2 sequestration now applies to the sum of all CO2 stores
  via an ``extra_functionality`` constraint.

* The myopic option can now be used together with different clustering for the
  generators and the network. The existing renewable capacities are split evenly
  among the regions in every country [`#144 <https://github.com/PyPSA/PyPSA-Eur-Sec/pull/144>`_].

* Add optional function to use ``geopy`` to locate entries of the Hotmaps
  database of industrial sites with missing location based on city and country,
  which reduces missing entries by half. It can be activated by setting
  ``industry: hotmaps_locate_missing: true``, takes a few minutes longer, and
  should only be used if spatial resolution is coarser than city level.


**Performance and Structure**

* Extended use of ``multiprocessing`` for much better performance
  (from up to 20 minutes to less than one minute).

* Handle most input files (or base directories) via ``snakemake.input``.

* Use of ``mock_snakemake`` from PyPSA-Eur.

* Update ``solve_network`` rule to match implementation in PyPSA-Eur by using
  ``n.ilopf()`` and remove outdated code using ``pyomo``.
  Allows the new setting to skip iterated impedance updates with ``solving:
  options: skip_iterations: true``.

* The component attributes that are to be overridden are now stored in the folder
  ``data/override_component_attrs`` analogous to ``pypsa/component_attrs``.
  This reduces verbosity and also allows circumventing the ``n.madd()`` hack
  for individual components with non-default attributes.
  This data is also tracked in the Snakefile.
  A function ``helper.override_component_attrs`` was added that loads this data
  and can pass the overridden component attributes into ``pypsa.Network()``.

* Add various parameters to ``config.default.yaml`` which were previously hardcoded inside the scripts
  (e.g. energy reference years, BEV settings, solar thermal collector models, geomap colours).

* Removed stale industry demand rules ``build_industrial_energy_demand_per_country``
  and ``build_industrial_demand``. These are superseded with more regionally resolved rules.

* Use simpler and shorter ``gdf.sjoin()`` function to allocate industrial sites
  from the Hotmaps database to onshore regions.
  This change also fixes a bug:
  The previous version allocated sites to the closest bus,
  but at country borders (where Voronoi cells are distorted by the borders),
  this had resulted in e.g. a Spanish site close to the French border
  being wrongly allocated to the French bus if the bus center was closer.

* Retrofitting rule is now only triggered if endogeneously optimised.

* Show progress in build rules with ``tqdm`` progress bars.

* Reduced verbosity of ``Snakefile`` through directory prefixes.

* Improve legibility of ``config.default.yaml`` and remove unused options.

* Use the country-specific time zone mappings from ``pytz`` rather than a manual mapping.

* A function ``add_carrier_buses()`` was added to the ``prepare_network`` rule to reduce code duplication.

* In the ``prepare_network`` rule the cost and potential adjustment was moved into an
  own function ``maybe_adjust_costs_and_potentials()``.

* Use ``matplotlibrc`` to set the default plotting style and backend.

* Added benchmark files for each rule.

* Consistent use of ``__main__`` block and further unspecific code cleaning.

* Updated data bundle and moved data bundle to zenodo.org (`10.5281/zenodo.5546517 <https://doi.org/10.5281/zenodo.5546517>`_).


**Bugfixes and Compatibility**

* Compatibility with ``atlite>=0.2``. Older versions of ``atlite`` will no longer work.

* Corrected calculation of "gas for industry" carbon capture efficiency.

* Implemented changes to ``n.snapshot_weightings`` in PyPSA v0.18.0.

* Compatibility with ``xarray`` version 0.19.

* New dependencies: ``tqdm``, ``atlite>=0.2.4``, ``pytz`` and ``geopy`` (optional).
  These are included in the environment specifications of PyPSA-Eur v0.4.0.

Many thanks to all who contributed to this release!


PyPSA-Eur-Sec 0.5.0 (21st May 2021)
===================================

This release includes improvements to the cost database for building retrofits, carbon budget management and wildcard settings, as well as an important bugfix for the emissions from land transport.

This release is known to work with `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ Version 0.3.0 and `Technology Data <https://github.com/PyPSA/technology-data>`_ Version 0.2.0.

Please note that the data bundle has also been updated.

New features and bugfixes:

* The cost database for retrofitting of the thermal envelope of buildings has been updated. Now, for calculating the space heat savings of a building, losses by thermal bridges and ventilation are included as well as heat gains (internal and by solar radiation). See the section :ref:`retro` for more details on the retrofitting module.
* For the myopic investment option, a carbon budget and a type of decay (exponential or beta) can be selected in the ``config.yaml`` file to distribute the budget across the ``planning_horizons``. For example, ``cb40ex0`` in the ``{sector_opts}`` wildcard will distribute a carbon budget of 40 GtCO2 following an exponential decay with initial growth rate 0.
* Added an option to alter the capital cost or maximum capacity of carriers by a factor via ``carrier+factor`` in the ``{sector_opts}`` wildcard. This can be useful for exploring uncertain cost parameters. Example: ``solar+c0.5`` reduces the ``capital_cost`` of solar to 50\% of original values. Similarly ``solar+p3`` multiplies the ``p_nom_max`` by 3.
* Rename the bus for European liquid hydrocarbons from ``Fischer-Tropsch`` to ``EU oil``, since it can be supplied not just with the Fischer-Tropsch process, but also with fossil oil.
* Bugfix: The new separation of land transport by carrier in Version 0.4.0 failed to account for the carbon dioxide emissions from internal combustion engines in land transport. This is now treated as a negative load on the atmospheric carbon dioxide bus, just like aviation emissions.
* Bugfix: Fix reading in of ``pypsa-eur/resources/powerplants.csv`` to PyPSA-Eur Version 0.3.0 (use column attribute name ``DateIn`` instead of old ``YearDecommissioned``).
* Bugfix: Make sure that ``Store`` components (battery and H2) are also removed from PyPSA-Eur, so they can be added later by PyPSA-Eur-Sec.

Thanks to Lisa Zeyen (KIT) for the retrofitting improvements and Marta Victoria (Aarhus University) for the carbon budget and wildcard management.

PyPSA-Eur-Sec 0.4.0 (11th December 2020)
=========================================

This release includes a more accurate nodal disaggregation of industry demand within each country, fixes to CHP and CCS representations, as well as changes to some configuration settings.

It has been released to coincide with `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ Version 0.3.0 and `Technology Data <https://github.com/PyPSA/technology-data>`_ Version 0.2.0, and is known to work with these releases.

New features:

* The `Hotmaps Industrial Database <https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`_ is used to disaggregate the industrial demand spatially to the nodes inside each country (previously it was distributed by population density).
* Electricity demand from industry is now separated from the regular electricity demand and distributed according to the industry demand. Only the remaining regular electricity demand for households and services is distributed according to GDP and population.
* A cost database for the retrofitting of the thermal envelope of residential and services buildings has been integrated, as well as endogenous optimisation of the level of retrofitting. This is described in the paper `Mitigating heat demand peaks in buildings in a highly renewable European energy system <https://arxiv.org/abs/2012.01831>`_. Retrofitting can be activated both exogenously and endogenously from the ``config.yaml``.
* The biomass and gas combined heat and power (CHP) parameters ``c_v`` and ``c_b`` were read in assuming they were extraction plants rather than back pressure plants. The data is now corrected in `Technology Data <https://github.com/PyPSA/technology-data>`_ Version 0.2.0 to the correct DEA back pressure assumptions and they are now implemented as single links with a fixed ratio of electricity to heat output (even as extraction plants, they were always sitting on the backpressure line in simulations, so there was no point in modelling the full heat-electricity feasibility polygon). The old assumptions underestimated the heat output.
* The Danish Energy Agency released `new assumptions for carbon capture <https://ens.dk/en/our-services/projections-and-models/technology-data/technology-data-industrial-process-heat-and>`_ in October 2020, which have now been incorporated in PyPSA-Eur-Sec, including direct air capture (DAC) and post-combustion capture on CHPs, cement kilns and other industrial facilities. The electricity and heat demand for DAC is modelled for each node (with heat coming from district heating), but currently the electricity and heat demand for industrial capture is not modelled very cleanly (for process heat, 10% of the energy is assumed to go to carbon capture) - a new issue will be opened on this.
* Land transport is separated by energy carrier (fossil, hydrogen fuel cell electric vehicle, and electric vehicle), but still needs to be separated into heavy and light vehicles (the data is there, just not the code yet).
* For assumptions that change with the investment year, there is a new time-dependent format in the ``config.yaml`` using a dictionary with keys for each year. Implemented examples include the CO2 budget, exogenous retrofitting share and land transport energy carrier; more parameters will be dynamised like this in future.
* Some assumptions have been moved out of the code and into the ``config.yaml``, including the carbon sequestration potential and cost, the heat pump sink temperature, reductions in demand for high value chemicals, and some BEV DSM parameters and transport efficiencies.
* Documentation on :doc:`supply_demand` options has been added.

Many thanks to Fraunhofer ISI for opening the hotmaps database and to Lisa Zeyen (KIT) for implementing the building retrofitting.


PyPSA-Eur-Sec 0.3.0 (27th September 2020)
=========================================

This releases focuses on improvements to industry demand and the generation of intermediate files for demand for basic materials. There are still inconsistencies with CCS and waste management that need to be improved.

It is known to work with PyPSA-Eur v0.1.0 (commit bb3477cd69), PyPSA v0.17.1 and technology-data v0.1.0. Please note that the data bundle has also been updated.


New features:

* In previous version of PyPSA-Eur-Sec the energy demand for industry was calculated directly for each location. Now, instead, the production of each material (steel, cement, aluminium) at each location is calculated as an intermediate data file, before the energy demand is calculated from it. This allows us in future to have competing industrial processes for supplying the same material demand.
* The script ``build_industrial_production_per_country_tomorrow.py`` determines the future industrial production of materials based on today's levels as well as assumed recycling and demand change measures.
* The energy demand for each industry sector and each location in 2015 is also calculated, so that it can be later incorporated in the pathway optimization.
* Ammonia production data is taken from the USGS and deducted from JRC-IDEES's "basic chemicals" so that it ammonia can be handled separately from the others (olefins, aromatics and chlorine).
* Solid biomass is no longer allowed to be used for process heat in cement and basic chemicals, since the wastes and residues cannot be guaranteed to reach the high temperatures required. Instead, solid biomass is used in the paper and pulp as well as food, beverages and tobacco industries, where required temperatures are lower (see `DOI:10.1002/er.3436 <https://doi.org/10.1002/er.3436>`_ and `DOI:10.1007/s12053-017-9571-y <https://doi.org/10.1007/s12053-017-9571-y>`_).
* National installable potentials for salt caverns are now applied.
* When electricity distribution grids are activated, new industry electricity demand, resistive heaters and micro-CHPs are now connected to the lower voltage levels.
* Gas distribution grid costs are included for gas boilers and micro-CHPs.
* Installable potentials for rooftop PV are included with an assumption of 1 kWp per person.
* Some intermediate files produced by scripts have been moved from the folder ``data`` to the folder ``resources``. Now ``data`` only includes input data, while ``resources`` only includes intermediate files necessary for building the network models. Please note that the data bundle has also been updated.
* Biomass potentials for different years and scenarios from the JRC are generated in an intermediate file, so that a selection can be made more explicitly by specifying the biomass types from the ``config.yaml``.


PyPSA-Eur-Sec 0.2.0 (21st August 2020)
======================================

This release introduces pathway optimization over many years (e.g. 2020, 2030, 2040, 2050) with myopic foresight, as well as outsourcing the technology assumptions to the `technology-data <https://github.com/PyPSA/technology-data>`_ repository.

It is known to work with PyPSA-Eur v0.1.0 (commit bb3477cd69), PyPSA v0.17.1 and technology-data v0.1.0.

New features:

* Option for pathway optimization with myopic foresight, based on the paper `Early decarbonisation of the European Energy system pays off (2020) <https://arxiv.org/abs/2004.11009>`_. Investments are optimized sequentially for multiple years (e.g. 2020, 2030, 2040, 2050) taking account of existing assets built in previous years and their lifetimes. The script uses data on the existing assets for electricity and building heating technologies, but there are no assumptions yet for existing transport and industry (if you include these, the model will greenfield them). There are also some `outstanding issues <https://github.com/PyPSA/pypsa-eur-sec/issues/19#issuecomment-678194802>`_ on e.g. the distribution of existing wind, solar and heating technologies within each country. To use myopic foresight, set ``foresight : 'myopic'`` in the ``config.yaml`` instead of the default ``foresight : 'overnight'``. An example configuration can be found in ``config.myopic.yaml``. More details on the implementation can be found in :doc:`myopic`.

* Technology assumptions (costs, efficiencies, etc.) are no longer stored in the repository. Instead, you have to install the `technology-data <https://github.com/PyPSA/technology-data>`_ database in a parallel directory. These assumptions are largely based on the `Danish Energy Agency Technology Data <https://ens.dk/en/our-services/projections-and-models/technology-data>`_. More details on the installation can be found in :doc:`installation`.

* Logs and benchmarks are now stored with the other model outputs in ``results/run-name/``.

* All buses now have a ``location`` attribute, e.g. bus ``DE0 3 urban central heat`` has a ``location`` of ``DE0 3``.

* All assets have a ``lifetime`` attribute (integer in years). For the myopic foresight, a ``build_year`` attribute is also stored.

* Costs for solar and onshore and offshore wind are recalculated by PyPSA-Eur-Sec based on the investment year, including the AC or DC connection costs for offshore wind.

Many thanks to Marta Victoria for implementing the myopic foresight, and Marta Victoria, Kun Zhu and Lisa Zeyen for developing the technology assumptions database.


PyPSA-Eur-Sec 0.1.0 (8th July 2020)
===================================

This is the first proper release of PyPSA-Eur-Sec, a model of the European energy system at the transmission network level that covers the full ENTSO-E area.

It is known to work with PyPSA-Eur v0.1.0 (commit bb3477cd69) and PyPSA v0.17.0.

We are making this release since in version 0.2.0 we will introduce changes to allow myopic investment planning that will require minor changes for users of the overnight investment planning.

PyPSA-Eur-Sec builds on the electricity generation and transmission
model `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ to add demand
and supply for the following sectors: transport, space and water
heating, biomass, industry and industrial feedstocks. This completes
the energy system and includes all greenhouse gas emitters except
waste management, agriculture, forestry and land use.

PyPSA-Eur-Sec was initially based on the model PyPSA-Eur-Sec-30 (Version 0.0.1 below) described
in the paper `Synergies of sector coupling and transmission
reinforcement in a cost-optimised, highly renewable European energy
system <https://arxiv.org/abs/1801.05290>`_ (2018) but it differs by
being based on the higher resolution electricity transmission model
`PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ rather than a
one-node-per-country model, and by including biomass, industry,
industrial feedstocks, aviation, shipping, better carbon management,
carbon capture and usage/sequestration, and gas networks.


PyPSA-Eur-Sec includes PyPSA-Eur as a
`snakemake <https://snakemake.readthedocs.io/en/stable/index.html>`_
`subworkflow <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-sub-workflows>`_. PyPSA-Eur-Sec
uses PyPSA-Eur to build the clustered transmission model along with
wind, solar PV and hydroelectricity potentials and time series. Then
PyPSA-Eur-Sec adds other conventional generators, storage units and
the additional sectors.




PyPSA-Eur-Sec 0.0.2 (4th September 2020)
========================================

This version, also called PyPSA-Eur-Sec-30-Path, built on
PyPSA-Eur-Sec 0.0.1 (also called PyPSA-Eur-Sec-30) to include myopic
pathway optimisation for the paper `Early decarbonisation of the
European energy system pays off <https://arxiv.org/abs/2004.11009>`_
(2020). The myopic pathway optimisation was then merged into the main
PyPSA-Eur-Sec codebase in Version 0.2.0 above.

This model has `its own github repository
<https://github.com/martavp/pypsa-eur-sec-30-path>`_ and is `archived
on Zenodo <https://zenodo.org/record/4014807>`_.



PyPSA-Eur-Sec 0.0.1 (12th January 2018)
========================================

This is the first published version of PyPSA-Eur-Sec, also called
PyPSA-Eur-Sec-30. It was first used in the research paper `Synergies of
sector coupling and transmission reinforcement in a cost-optimised,
highly renewable European energy system
<https://arxiv.org/abs/1801.05290>`_ (2018). The model covers 30
European countries with one node per country. It includes demand and
supply for electricity, space and water heating in buildings, and land
transport.

It is `archived on Zenodo <https://zenodo.org/record/1146666>`_.



Release Process
===============

* Finalise release notes at ``doc/release_notes.rst``.

* Update version number in ``doc/conf.py`` and ``*config.*.yaml``.

* Make a ``git commit``.

* Tag a release by running ``git tag v0.x.x``, ``git push``, ``git push --tags``. Include release notes in the tag message.

* Make a `GitHub release <https://github.com/PyPSA/pypsa-eur-sec/releases>`_, which automatically triggers archiving by `zenodo <https://doi.org/10.5281/zenodo.3938042>`_.

* Send announcement on the `PyPSA mailing list <https://groups.google.com/forum/#!forum/pypsa>`_.

To make a new release of the data bundle, make an archive of the files in ``data`` which are not already included in the git repository:

.. code:: bash

    data % tar pczf pypsa-eur-sec-data-bundle.tar.gz eea/UNFCCC_v23.csv switzerland-sfoe biomass eurostat-energy_balances-* jrc-idees-2015 emobility WindWaveWEC_GLTB.xlsx myb1-2017-nitro.xls Industrial_Database.csv retro/tabula-calculator-calcsetbuilding.csv nuts/NUTS_RG_10M_2013_4326_LEVL_2.geojson h2_salt_caverns_GWh_per_sqkm.geojson
