##########################################
Release Notes
##########################################

Future release
===================

* For the myopic investment option, a carbon budget and a type of decay (exponential or beta) can be selected in the ``config.yaml`` file to distribute the budget across the ``planning_horizons``. For example, ``cb40ex0`` in the ``{sector_opts}`` wildcard will distribute a carbon budget of 40 GtCO2 following an exponential decay with initial growth rate 0.
* Added an option to alter the capital cost or maximum capacity of carriers by a factor via ``carrier+factor`` in the ``{sector_opts}`` wildcard. This can be useful for exploring uncertain cost parameters. Example: ``solar+c0.5`` reduces the ``capital_cost`` of solar to 50\% of original values. Similarly ``solar+p3`` multiplies the ``p_nom_max`` by 3.
* Rename the bus for European liquid hydrocarbons from ``Fischer-Tropsch`` to ``EU oil``, since it can be supplied not just with the Fischer-Tropsch process, but also with fossil oil.
* Bugfix: Fix reading in of ``pypsa-eur/resources/powerplants.csv`` to PyPSA-Eur Version 0.3.0 (use column attribute name ``DateIn`` instead of old ``YearDecommissioned``).
* Bugfix: Make sure that ``Store`` components (battery and H2) are also removed from PyPSA-Eur, so they can be added later by PyPSA-Eur-Sec.


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

This is the first release of PyPSA-Eur-Sec, a model of the European energy system at the transmission network level that covers the full ENTSO-E area.

It is known to work with PyPSA-Eur v0.1.0 (commit bb3477cd69) and PyPSA v0.17.0.

We are making this release since in version 0.2.0 we will introduce changes to allow myopic investment planning that will require minor changes for users of the overnight investment planning.

PyPSA-Eur-Sec builds on the electricity generation and transmission
model `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ to add demand
and supply for the following sectors: transport, space and water
heating, biomass, industry and industrial feedstocks. This completes
the energy system and includes all greenhouse gas emitters except
waste management, agriculture, forestry and land use.

PyPSA-Eur-Sec was initially based on the model PyPSA-Eur-Sec-30 described
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

    data % tar pczf pypsa-eur-sec-data-bundle-YYMMDD.tar.gz eea/UNFCCC_v23.csv switzerland-sfoe biomass eurostat-energy_balances-* jrc-idees-2015 emobility urban_percent.csv timezone_mappings.csv heat_load_profile_DK_AdamJensen.csv WindWaveWEC_GLTB.xlsx myb1-2017-nitro.xls Industrial_Database.csv
