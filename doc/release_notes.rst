..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Release Notes
##########################################

Upcoming Release
================

* Bugfix: Correct technology keys for the electricity production plotting to work out the box.

* New configuration option ``everywhere_powerplants`` to build conventional powerplants everywhere, irrespective of existing powerplants locations, in the network (https://github.com/PyPSA/pypsa-eur/pull/850).

* Remove option for wave energy as technology data is not maintained.

* Define global constraint for CO2 emissions on the final state of charge of the
  CO2 atmosphere store. This gives a more sparse constraint that should improve
  the performance of the solving process.

* Bugfix: Assure entering of code block which corrects Norwegian heat demand.

* Add warning when BEV availability weekly profile has negative values in `build_transport_demand`.

* Stacktrace of uncaught exceptions should now be correctly included inside log files (via `configure_logging(..)`).

* Cluster residential and services heat buses by default. Can be disabled with ``cluster_heat_buses: false``.

* Bugfix: Do not reduce district heat share when building population-weighted
  energy statistics. Previously the district heating share was being multiplied
  by the population weighting, reducing the DH share with multiple nodes.

* Move building of daily heat profile to its own rule
  :mod:`build_hourly_heat_demand` from :mod:`prepare_sector_network`.

* In :mod:`build_energy_totals`, district heating shares are now reported in a
  separate file.

* Move calculation of district heating share to its own rule
  :mod:`build_district_heat_share`.

* Move building of distribution of existing heating to own rule
  :mod:`build_existing_heating_distribution`. This makes the distribution of
  existing heating to urban/rural, residential/services and spatially more
  transparent.

* Bugfix: Correctly read out number of solver threads from configuration file.

* Air-sourced heat pumps can now also be built in rural areas. Previously, only
  ground-sourced heat pumps were considered for this category.

* Bugfix: Correctly read out number of solver threads from configuration file.

* Add support for the linopy ``io_api`` option; set to ``"direct"`` to increase model reading and writing performance for the highs and gurobi solvers.

* Add the option to customise map projection in plotting config.

* The rule ``plot_network`` has been split into separate rules for plotting
  electricity, hydrogen and gas networks.


PyPSA-Eur 0.9.0 (5th January 2024)
==================================

**New Features**

* Add option to specify losses for bidirectional links, e.g. pipelines or HVDC
  links, in configuration file under ``sector: transmission_efficiency:``. Users
  can specify static or length-dependent values as well as a length-dependent
  electricity demand for compression, which is implemented as a multi-link to
  the local electricity buses. The bidirectional links will then be split into
  two unidirectional links with linked capacities (https://github.com/PyPSA/pypsa-eur/pull/739).

* Merged option to extend geographical scope to Ukraine and Moldova. These
  countries are excluded by default and is currently constrained to power-sector
  only parts of the workflow. A special config file
  `config/config.entsoe-all.yaml` was added as an example to run the workflow
  with all ENTSO-E member countries (including observer members like Ukraine and
  Moldova). Moldova can currently only be included in conjunction with Ukraine
  due to the absence of demand data. The Crimean power system is manually
  reconnected to the main Ukrainian grid with the configuration option
  `reconnect_crimea` (https://github.com/PyPSA/pypsa-eur/pull/321).

* New experimental support for multi-decade optimisation with perfect foresight
  (``foresight: perfect``). Maximum growth rates for carriers, global carbon
  budget constraints and emission constraints for particular investment periods.

* Add option to reference an additional source file where users can specify
  custom ``extra_functionality`` constraints in the configuration file. The
  default setting points to an empty hull at
  ``data/custom_extra_functionality.py`` (https://github.com/PyPSA/pypsa-eur/pull/824).

* Add locations, capacities and costs of existing gas storage using Global
  Energy Monitor's `Europe Gas Tracker
  <https://globalenergymonitor.org/projects/europe-gas-tracker>`_
  (https://github.com/PyPSA/pypsa-eur/pull/835).

* Add option to use `LUISA Base Map
  <https://publications.jrc.ec.europa.eu/repository/handle/JRC124621>`_ 50m land
  coverage dataset for land eligibility analysis in
  :mod:`build_renewable_profiles`. Settings are analogous to the CORINE dataset
  but with the key ``luisa:`` in the configuration file. To leverage the
  dataset's full advantages, set the excluder resolution to 50m
  (``excluder_resolution: 50``). For land category codes, see `Annex 1 of the
  technical documentation
  <https://publications.jrc.ec.europa.eu/repository/bitstream/JRC124621/technical_report_luisa_basemap_2018_v7_final.pdf>`_
  (https://github.com/PyPSA/pypsa-eur/pull/842).

* Add option to capture CO2 contained in biogas when upgrading (``sector:
  biogas_to_gas_cc``) (https://github.com/PyPSA/pypsa-eur/pull/615).

* If load shedding is activated, it is now applied to all carriers, not only
  electricity (https://github.com/PyPSA/pypsa-eur/pull/784).

* Add option for heat vents in district heating (``sector:
  central_heat_vent:``). The combination of must-run conditions for some
  power-to-X processes, waste heat usage enabled and decreasing heating demand,
  can lead to infeasibilities in pathway optimisation for some investment
  periods since larger Fischer-Tropsch capacities are needed in early years but
  the waste heat exceeds the heat demand in later investment periods.
  (https://github.com/PyPSA/pypsa-eur/pull/791).

* Allow possibility to go from copperplated to regionally resolved methanol and
  oil demand with switches ``sector: regional_methanol_demand: true`` and
  ``sector: regional_oil_demand: true``. This allows nodal/regional CO2
  constraints to be applied (https://github.com/PyPSA/pypsa-eur/pull/827).

* Allow retrofitting of existing gas boilers to hydrogen boilers in pathway
  optimisation.

* Add option to add time-varying CO2 emission prices (electricity-only, ``costs:
  emission_prices: co2_monthly_prices: true``). This is linked to the new
  ``{opts}`` wildcard option ``Ept``.

* Network clustering can now consider efficiency classes when aggregating
  carriers. The option ``clustering: consider_efficiency_classes:`` aggregates
  each carriers into the top 10-quantile (high), the bottom 90-quantile (low),
  and everything in between (medium).

* Added option ``conventional: dynamic_fuel_price:`` to consider the monthly
  fluctuating fuel prices for conventional generators. Refer to the CSV file
  ``data/validation/monthly_fuel_price.csv``.

* For hydro-electricity, add switches ``flatten_dispatch`` to consider an upper
  limit for the hydro dispatch. The limit is given by the average capacity
  factor plus the buffer given in  ``flatten_dispatch_buffer``.

* Extend options for waste heat usage from Haber-Bosch, methanolisation and
  methanation (https://github.com/PyPSA/pypsa-eur/pull/834).

* Add new ``sector_opts`` wildcard option "nowasteheat" to disable all waste
  heat usage (https://github.com/PyPSA/pypsa-eur/pull/834).

* Add new rule ``retrieve_irena`` to automatically retrieve up-to-date values
  for existing renewables capacities (https://github.com/PyPSA/pypsa-eur/pull/756).

* Print Irreducible Infeasible Subset (IIS) if model is infeasible. Only for
  solvers with IIS support (https://github.com/PyPSA/pypsa-eur/pull/841).

* More wildcard options now have a corresponding config entry. If the wildcard
  is given, then its value is used. If the wildcard is not given but the options
  in config are enabled, then the value from config is used. If neither is
  given, the options are skipped (https://github.com/PyPSA/pypsa-eur/pull/827).

* Validate downloads from Zenodo using MD5 checksums. This identifies corrupted
  or incomplete downloads (https://github.com/PyPSA/pypsa-eur/pull/821).

* Add rule ``sync`` to synchronise with a remote machine using the ``rsync``
  library. Configuration settings are found under ``remote:``.

**Breaking Changes**

* Remove all negative loads on the ``co2 atmosphere`` bus representing emissions
  for e.g. fixed fossil demands for transport oil. Instead these are handled
  more transparently with a fixed transport oil demand and a link taking care of
  the emissions to the ``co2 atmosphere`` bus. This is also a preparation for
  endogenous transport optimisation, where demand will be subject to
  optimisation (e.g. fuel switching in the transport sector)
  (https://github.com/PyPSA/pypsa-eur/pull/827).

* Process emissions from steam crackers (i.e. naphtha processing for HVC) are
  now piped from the consumption link to the process emissions bus where the
  model can decide about carbon capture. Previously the process emissions for
  naphtha were a fixed load (https://github.com/PyPSA/pypsa-eur/pull/827).

* Distinguish between stored and sequestered CO2. Stored CO2 is stored
  overground in tanks and can be used for CCU (e.g. methanolisation).
  Sequestered CO2 is stored underground and can no longer be used for CCU. This
  distinction is made because storage in tanks is more expensive than
  underground storage. The link that connects stored and sequestered CO2 is
  unidirectional (https://github.com/PyPSA/pypsa-eur/pull/844).

* Files extracted from sector-coupled data bundle have been moved from ``data/``
  to ``data/sector-bundle``.

* Split configuration to enable SMR and SMR CC (``sector: smr:`` and ``sector:
  smr_cc:``) (https://github.com/PyPSA/pypsa-eur/pull/757).

* Add separate option to add resistive heaters to the technology choices
  (``sector: resistive_heaters:``). Previously they were always added when
  boilers were added (https://github.com/PyPSA/pypsa-eur/pull/808).

* Remove HELMETH option (``sector: helmeth:``).

* Remove "conservative" renewable potentials estimation option
  (https://github.com/PyPSA/pypsa-eur/pull/838).

* With this release we stop posting updates to the network pre-builts.

**Changes**

* Updated Global Energy Monitor LNG terminal data to March 2023 version
  (https://github.com/PyPSA/pypsa-eur/pull/707).

* For industry distribution, use EPRTR as fallback if ETS data is not available
  (https://github.com/PyPSA/pypsa-eur/pull/721).

* It is now possible to specify years for biomass potentials which do not exist
  in the JRC-ENSPRESO database, e.g. 2037. These are linearly interpolated
  (https://github.com/PyPSA/pypsa-eur/pull/744).

* In pathway mode, the biomass potential is linked to the investment year
  (https://github.com/PyPSA/pypsa-eur/pull/744).

* Increase allowed deployment density of solar to 5.1 MW/sqkm by default.

* Default to full electrification of land transport by 2050.

* Provide exogenous transition settings in 5-year steps.

* Default to approximating transmission losses in HVAC lines
  (``transmission_losses: 2``).

* Use electrolysis waste heat by default.

* Set minimum part loads for PtX processes to 30% for methanolisation and
  methanation, and to 70% for Fischer-Tropsch synthesis.

* Add VOM as marginal cost to PtX processes
  (https://github.com/PyPSA/pypsa-eur/pull/830).

* Add pelletizing costs for biomass boilers (https://github.com/PyPSA/pypsa-eur/pull/833).

* Update default offshore wind turbine model to "NREL Reference 2020 ATB 5.5 MW"
  (https://github.com/PyPSA/pypsa-eur/pull/832).

* Switch to using hydrogen and electricity inputs for Haber-Bosch from
  https://github.com/PyPSA/technology-data (https://github.com/PyPSA/pypsa-eur/pull/831).

* The configuration setting for country focus weights when clustering the
  network has been moved from ``focus_weights:`` to ``clustering:
  focus_weights:``. Backwards compatibility to old config files is maintained
  (https://github.com/PyPSA/pypsa-eur/pull/794).

* The ``mock_snakemake`` function can now be used with a Snakefile from a
  different directory using the new ``root_dir`` argument
  (https://github.com/PyPSA/pypsa-eur/pull/771).

* Rule ``purge`` now initiates a dialog to confirm if purge is desired
  (https://github.com/PyPSA/pypsa-eur/pull/745).

* Files downloaded from zenodo are now write-protected to prevent accidental
  re-download (https://github.com/PyPSA/pypsa-eur/pull/730).

* Performance improvements for rule ``build_ship_raster``
  (https://github.com/PyPSA/pypsa-eur/pull/845).

* Improve time logging in :mod:`build_renewable_profiles`
  (https://github.com/PyPSA/pypsa-eur/pull/837).

* In myopic pathway optimisation, disable power grid expansion if line volume
  already hit (https://github.com/PyPSA/pypsa-eur/pull/840).

* JRC-ENSPRESO data is now downloaded from a Zenodo mirror because the link was
  unreliable (https://github.com/PyPSA/pypsa-eur/pull/801).

* Add focus weights option for clustering to documentation
  (https://github.com/PyPSA/pypsa-eur/pull/781).

* Add proxy for biomass transport costs if no explicit biomass transport network
  is considered (https://github.com/PyPSA/pypsa-eur/pull/711).

**Bugs and Compatibility**

* The minimum PyPSA version is now 0.26.1.

* Update to ``tsam>=0.2.3`` for performance improvents in temporal clustering.

* Pin ``snakemake`` version to below 8.0.0, as the new version is not yet
  supported. The next release will switch to the requirement ``snakemake>=8``.

* Bugfix: Add coke and coal demand for integrated steelworks
  (https://github.com/PyPSA/pypsa-eur/pull/718).

* Bugfix: Make :mod:`build_renewable_profiles` consider subsets of cutout time
  scope (https://github.com/PyPSA/pypsa-eur/pull/709).

* Bugfix: In :mod:`simplify network`, remove 'underground' column to avoid
  consense error (https://github.com/PyPSA/pypsa-eur/pull/714).

* Bugfix: Fix in :mod:`add_existing_baseyear` to account for the case when there
  is no rural heating demand for some nodes in network
  (https://github.com/PyPSA/pypsa-eur/pull/706).

* Bugfix: The unit of the capital cost of Haber-Bosch plants was corrected
  (https://github.com/PyPSA/pypsa-eur/pull/829).

* The minimum capacity for renewable generators when using the myopic option has
  been fixed (https://github.com/PyPSA/pypsa-eur/pull/728).

* Compatibility for running with single node and single country
  (https://github.com/PyPSA/pypsa-eur/pull/839).

* A bug preventing the addition of custom powerplants specified in
  ``data/custom_powerplants.csv`` was fixed.
  (https://github.com/PyPSA/pypsa-eur/pull/732)

* Fix nodal fraction in :mod:`add_existing_year` when using distributed
  generators (https://github.com/PyPSA/pypsa-eur/pull/798).

* Bugfix: District heating without progress caused division by zero
  (https://github.com/PyPSA/pypsa-eur/pull/796).

* Bugfix: Drop duplicates in :mod:`build_industrial_distribution_keys`, which
  can occur through the geopandas ``.sjoin()`` function if a point is located on
  a border (https://github.com/PyPSA/pypsa-eur/pull/726).

* For network clustering fall back to ``ipopt`` when ``highs`` is designated
  solver (https://github.com/PyPSA/pypsa-eur/pull/795).

* Fix typo in buses definition for oil boilers in ``add_industry`` in
  :mod:`prepare_sector_network` (https://github.com/PyPSA/pypsa-eur/pull/812).

* Resolve code issues for endogenous building retrofitting. Select correct
  sector names, address deprecations, distinguish between district heating,
  decentral heating in urban areas or rural areas for floor area calculations
  (https://github.com/PyPSA/pypsa-eur/pull/808).

* Addressed various deprecations.


PyPSA-Eur 0.8.1 (27th July 2023)
================================

**New Features**

* Add option to consider dynamic line rating based on wind speeds and
  temperature according to `Glaum and Hofmann (2022)
  <https://arxiv.org/abs/2208.04716>`_. See configuration section ``lines:
  dynamic_line_rating:`` for more details. (https://github.com/PyPSA/pypsa-eur/pull/675)

* Add option to include a piecewise linear approximation of transmission losses,
  e.g. by setting ``solving: options: transmission_losses: 2`` for an
  approximation with two tangents. (https://github.com/PyPSA/pypsa-eur/pull/664)

* Add plain hydrogen turbine as additional re-electrification option besides
  hydrogen fuel cell. Add switches for both re-electrification options under
  ``sector: hydrogen_turbine:`` and ``sector: hydrogen_fuel_cell:``.
  (https://github.com/PyPSA/pypsa-eur/pull/647)

* Added configuration option ``lines: max_extension:`` and ``links:
  max_extension:``` to control the maximum capacity addition per line or link in
  MW. (https://github.com/PyPSA/pypsa-eur/pull/665)

* A ``param:`` section in the snakemake rule definitions was added to track
  changed settings in ``config.yaml``. The goal is to automatically re-execute
  rules where parameters have changed. See `Non-file parameters for rules
  <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#non-file-parameters-for-rules>`_
  in the snakemake documentation. (https://github.com/PyPSA/pypsa-eur/pull/663)

* A new function named ``sanitize_carrier`` ensures that all unique carrier
  names are present in the network's carriers attribute, and adds nice names and
  colors for each carrier according to the provided configuration dictionary.
  (https://github.com/PyPSA/pypsa-eur/pull/653,
  https://github.com/PyPSA/pypsa-eur/pull/690)

* The configuration settings have been documented in more detail.
  (https://github.com/PyPSA/pypsa-eur/pull/685)

**Breaking Changes**

* The configuration files are now located in the ``config`` directory. This
  includes the ``config.default.yaml``, ``config.yaml`` as well as the test
  configuration files which are now located in the ``config/test`` directory.
  Config files that are still in the root directory will be ignored.
  (https://github.com/PyPSA/pypsa-eur/pull/640)

* Renamed script and rule name from ``build_load_data`` to
  ``build_electricity_demand`` and ``retrieve_load_data`` to
  ``retrieve_electricity_demand``. (https://github.com/PyPSA/pypsa-eur/pull/642,
  https://github.com/PyPSA/pypsa-eur/pull/652)

* Updated to new spatial clustering module introduced in PyPSA v0.25.
  (https://github.com/PyPSA/pypsa-eur/pull/696)

**Changes**

* Handling networks with links with multiple inputs/outputs no longer requires
  to override component attributes.
  (https://github.com/PyPSA/pypsa-eur/pull/695)

* Added configuration option ``enable: retrieve:`` to control whether data
  retrieval rules from snakemake are enabled or not. Th default setting ``auto``
  will automatically detect and enable/disable the rules based on internet
  connectivity. (https://github.com/PyPSA/pypsa-eur/pull/694)

* Update to ``technology-data`` v0.6.0.
  (https://github.com/PyPSA/pypsa-eur/pull/704)

* Handle data bundle extraction paths via ``snakemake.output``.

* Additional technologies are added to ``tech_color`` in the configuration files
  to include previously unlisted carriers.

* Doc: Added note that Windows is only tested in CI with WSL.
  (https://github.com/PyPSA/pypsa-eur/issues/697)

* Doc: Add support section. (https://github.com/PyPSA/pypsa-eur/pull/656)

* Open ``rasterio`` files with ``rioxarray``.
  (https://github.com/PyPSA/pypsa-eur/pull/474)

* Migrate CI to ``micromamba``. (https://github.com/PyPSA/pypsa-eur/pull/700)

**Bugs and Compatibility**

* The new minimum PyPSA version is v0.25.1.

* Removed ``vresutils`` dependency.
  (https://github.com/PyPSA/pypsa-eur/pull/662)

* Adapt to new ``powerplantmatching`` version.
  (https://github.com/PyPSA/pypsa-eur/pull/687,
  https://github.com/PyPSA/pypsa-eur/pull/701)

* Bugfix: Correct typo in the CPLEX solver configuration in
  ``config.default.yaml``. (https://github.com/PyPSA/pypsa-eur/pull/630)

* Bugfix: Error in ``add_electricity`` where carriers were added multiple times
  to the network, resulting in a non-unique carriers error.

* Bugfix of optional reserve constraint.
  (https://github.com/PyPSA/pypsa-eur/pull/645)

* Fix broken equity constraints logic.
  (https://github.com/PyPSA/pypsa-eur/pull/679)

* Fix addition of load shedding generators.
  (https://github.com/PyPSA/pypsa-eur/pull/649)

* Fix automatic building of documentation on readthedocs.org.
  (https://github.com/PyPSA/pypsa-eur/pull/658)

* Bugfix: Update network clustering to avoid adding deleted links in clustered
  network. (https://github.com/PyPSA/pypsa-eur/pull/678)

* Address ``geopandas`` deprecations.
  (https://github.com/PyPSA/pypsa-eur/pull/678)

* Fix bug with underground hydrogen storage creation, where for some small model
  regions no cavern storage is available.
  (https://github.com/PyPSA/pypsa-eur/pull/672)


* Addressed deprecation warnings for ``pandas=2.0``. ``pandas=2.0`` is now minimum requirement.

PyPSA-Eur 0.8.0 (18th March 2023)
=================================

.. note::
  This is the first release of PyPSA-Eur which incorporates its sector-coupled extension PyPSA-Eur-Sec (v0.7.0).
  PyPSA-Eur can now directly be used for high-resolution energy system modelling with sector-coupling
  including industry, transport, buildings, biomass, and detailed carbon management. The PyPSA-Eur-Sec repository is now deprecated.

* The :mod:`solve_network` script now uses the ``linopy`` backend of PyPSA and is applied for both electricity-only and sector-coupled models. This
  requires an adjustment of custom ``extra_functionality``.
  See the `migration guide <https://pypsa.readthedocs.io/en/latest/examples/optimization-with-linopy-migrate-extra-functionalities.html>`_ in the PyPSA documentation.

* The configuration file ``config.default.yaml`` now also includes settings for
  sector-coupled models, which will be ignored when the user runs
  electricity-only studies. Common settings have been aligned.

* Unified handling of scenario runs. Users can name their scenarios in ``run:
  name:``, which will encapsulate results in a correspondingly named folder
  under ``results``. Additionally, users can select to encapsulate the ``resources`` folder
  in the same way, through the setting ``run: shared_resources:``.

* The solver configurations in ``config.default.yaml`` are now modularized. To
  change the set of solver options, change to value in ``solving: solver:
  options:`` to one of the keys in ``solving: solver_options:``.

* The ``Snakefile`` has been modularised. Rules are now organised in the
  ``rules`` directory.

* Unified wildcard for transmission line expansion from ``{lv}`` and ``{ll}`` to
  ``{ll}``.

* Renamed collection rules to distinguish between sector-coupled and
  electricity-only runs: ``cluster_networks``, ``extra_components_networks``,
  ``prepare_elec_networks``, ``prepare_sector_networks``,
  ``solve_elec_networks``, ``solve_sector_networks``, ``plot_networks``,
  ``all``.

* Some rules with a small computational footprint have been declared as ``localrules``.

* Added new utility rules ``purge`` for clearing workflow outputs from the
  directory, ``doc`` to build the documentation, and ``dag`` to create a
  workflow graph.

* The workflow can now be used with the ``snakemake --use-conda`` directive. In
  this way, Snakemake can automatically handle the installation of dependencies.

* Data retrieval rules now retry download twice in case of connection problems.

* The cutouts are now marked as ``protected()`` in the workflow to avoid
  accidental recomputation.

* The files contained in ``data/bundle`` are now marked as ``ancient()`` as they
  are not expected to be altered by workflow changes.

* Preparation scripts for sector-coupled models have been improved to only run
  for the subset of selected countries rather than all European countries.

* Added largely automated country code conversion using ``country_converter``..

* Test coverage extended to an electricity-only run and sector-coupled runs for
  overnight and myopic foresight scenarios for Ubuntu, MacOS and Windows.

* Apply ``black`` and ``snakefmt`` code formatting.

* Implemented REUSE compatibility for merged code.

* Merged documentations of PyPSA-Eur and PyPSA-Eur-Sec.

* Added a tutorial for running sector-coupled models to the documentation
  (:ref:`tutorial_sector`).

* Deleted ``config.tutorial.yaml``, which is superseded by
  ``test/config.electricity.yaml``.

* The ``mock_snakemake`` function now also takes configuration files as inputs.

* The helper scripts ``helper.py`` and ``_helpers.py`` have been merged into
  ``_helpers.py``.

* The unused rule ``plot_p_nom_max`` has been removed.

* The rule ``solve_network`` from PyPSA-Eur-Sec was renamed to
  ``solve_sector_network``.

* The plotting scripts from PyPSA-Eur (electricity-only) have been removed and
  are superseded by those from PyPSA-Eur-Sec (sector-coupled).

PyPSA-Eur Releases (pre-merge)
==============================

PyPSA-Eur 0.7.0 (16th February 2023)
------------------------------------


**New Features**

* Carriers of generators can now be excluded from aggregation in clustering
  network and simplify network (see ``exclude_carriers``).

* Added control for removing stubs in  :mod:`simplify_network` with options
  ``remove_stubs`` and ``remove_stubs_across_countries``.

* Add control for showing a progressbar in ``atlite`` processes
  (``show_progress``). Disabling the progressbar saves a lot of time.

* Added control for resolution of land eligibility analysis (see
  ``excluder_resolution``).


**Breaking Changes**

* The config entry ``snapshots: closed:`` was renamed to ``snapshots:
  inclusive:`` to address the upstream deprecation with ``pandas=1.4``. The
  previous setting ``None`` is no longer supported and replaced by ``both``, see
  the `pandas documentation
  <https://pandas.pydata.org/docs/reference/api/pandas.date_range.html>`_.
  Minimum version is now ``pandas>=1.4``.

* The configuration setting ``summary_dir`` was removed.


**Changes**

* Configuration defaults to new ``technology-data`` version 0.5.0.

* Fixed CRS warnings when projection of datasets was not specified.

* Cleaned shape unary unions.

* Increased resource requirements for some rules.

* Updated documentation.

* The documentation now uses the ``sphinx_book_theme``.


**Bugs and Compatibility**


* Bugfix: Corrected extent of natural protection areas in :mod:`build_natura_raster`.

* Bugfix: Use correct load variables for formulating reserve constraints.

* Bugfix: Use all available energy-to-power ratios for hydropower plants.

* Bugfix: The most recent processing of the ``entsoegridkit`` extract required
  further manual corrections. Also, the connection points of TYNDP links were
  corrected.

* Bugfix: Handle absence of hydropower inflow in ``EQ`` constraint.

* Compatibility with ``pyomo>=6.4.3`` in :mod:`cluster_network`.

* Upgrade to ``shapely>=2``.

* Updated version of CI cache action to version 3.
*
* Updated dependency constraints in ``environment.yaml``.

* Address various deprecation warnings.



PyPSA-Eur 0.6.1 (20th September 2022)
-------------------------------------

* Individual commits are now tested against pre-commit hooks. This includes
  black style formatting, sorting of package imports, Snakefile formatting and
  others. Installation instructions can for the pre-commit can be found `here
  <https://pre-commit.com/>`_.

* Pre-commit CI is now part of the repository's CI.

* The software now supports running the workflow with different settings within
  the same directory. A new config section ``run`` was created that specifies
  under which scenario ``name`` the created resources, networks and results
  should be stored. If ``name`` is not specified, the workflow uses the default
  paths. The entry ``shared_cutouts`` specifies whether the run should use
  cutouts from the default root directory or use run-specific cutouts.

* The heuristic distribution of today's renewable capacity installations is now
  enabled by default.

* The marginal costs of conventional generators are now taking the plant-specific
  efficiency into account where available.

PyPSA-Eur 0.6.0 (10th September 2022)
-------------------------------------

* Functionality to consider shipping routes when calculating the available area
  for offshore technologies were added. Data for the shipping density comes from
  the `Global Shipping Traffic Density dataset
  <https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`_.

* When transforming all transmission lines to a unified voltage level of 380kV,
  the workflow now preserves the transmission capacity rather than electrical
  impedance and reactance.

* Memory resources are now specified for all rules.

* Filtering of power plant data was adjusted to new versions of
  ``powerplantmatching``.

* The resolution of land exclusion calculation is now a configurable option. See
  setting ``excluder_resolution``.


PyPSA-Eur 0.5.0 (27th July 2022)
--------------------------------

**New Features**

* New network topology extracted from the ENTSO-E interactive map.
* Added existing renewable capacities for all countries based on IRENA
  statistics (IRENASTAT) using new ``powerplantmatching`` version:
* The corresponding ``config`` entries changed from ``estimate_renewable_capacities_from_capacity_stats`` to ``estimate_renewable_capacities``.
* The estimation is endabled by setting the subkey ``enable`` to ``True``.
* Configuration of reference year for capacities can be configured (default: ``2020``)
* The list of renewables provided by the OPSD database can be used as a basis, using the tag ``from_opsd: True``. This adds the renewables from the database and fills up the missing capacities with the heuristic distribution.
* Uniform expansion limit of renewable build-up based on existing capacities
  can be configured using ``expansion_limit`` option (default: ``false``;
  limited to determined renewable potentials)
* Distribution of country-level capacities proportional to maximum annual
  energy yield for each bus region
* The config key ``renewable_capacities_from_OPSD`` is deprecated and was moved
  under the section, ``estimate_renewable_capacities``. To enable it, set
  ``from_opsd`` to ``True``.

* Add operational reserve margin constraint analogous to `GenX implementation
  <https://genxproject.github.io/GenX/dev/core/#Reserves>`_. Can be activated
  with config setting ``electricity: operational_reserve:``.

* Implement country-specific  Energy Availability Factors (EAFs) for nuclear
  power plants based on IAEA 2018-2020 reported country averages. These are
  specified ``data/nuclear_p_max_pu.csv`` and translate to static ``p_max_pu``
  values.

* Add function to add global constraint on use of gas in :mod:`prepare_network`.
  This can be activated by including the keyword ``CH4L`` in the ``{opts}``
  wildcard which enforces the limit set in ``electricity: gaslimit:`` given in
  MWh thermal. Alternatively, it is possible to append a number in the ``{opts}``
  wildcard, e.g. ``CH4L200`` which limits the gas use to 200 TWh thermal.

* Add option to alter marginal costs of a carrier through ``{opts}`` wildcard:
  ``<carrier>+m<factor>``, e.g. ``gas+m2.5``, will multiply the default marginal
  cost for gas by factor 2.5.

* Hierarchical clustering was introduced. Distance metric is calculated from
  renewable potentials on hourly (feature entry ends with ``-time``) or annual
  (feature entry in config end with ``-cap``) values.

* Greedy modularity clustering was introduced. Distance metric is based on electrical distance taking into account the impedance of all transmission lines of the network.

* Techno-economic parameters of technologies (e.g. costs and efficiencies) will
  now be retrieved from a separate repository `PyPSA/technology-data
  <https://github.com/pypsa/technology-data>`_ that collects assumptions from a
  variety of sources. It is activated by default with ``enable:
  retrieve_cost_data: true`` and controlled with ``costs: year:`` and ``costs:
  version:``. The location of this data changed from ``data/costs.csv`` to
  ``resources/costs.csv`` [`#184
  <https://github.com/PyPSA/pypsa-eur/pull/184>`_].

* A new section ``conventional`` was added to the config file. This section
  contains configurations for conventional carriers.

* Add configuration option to implement arbitrary generator attributes for
  conventional generation technologies.

* Add option to set CO2 emission prices through ``{opts}`` wildcard: ``Ep<number>``,
  e.g. ``Ep180``, will set the EUR/tCO2 price.

**Changes**

* Add an efficiency factor of 88.55% to offshore wind capacity factors as a
  proxy for wake losses. More rigorous modelling is `planned
  <https://github.com/PyPSA/pypsa-eur/issues/153>`_ [`#277
  <https://github.com/PyPSA/pypsa-eur/pull/277>`_].

* Following discussion in `#285
  <https://github.com/PyPSA/pypsa-eur/issues/285>`_ we have disabled the
  correction factor for solar PV capacity factors by default while satellite
  data is used. A correction factor of 0.854337 is recommended if reanalysis
  data like ERA5 is used.

* The default deployment density of AC- and DC-connected offshore wind capacity
  is reduced from 3 MW/sqkm to a more conservative estimate of 2 MW/sqkm [`#280
  <https://github.com/PyPSA/pypsa-eur/pull/280>`_].

* The inclusion of renewable carriers is now specified in the config entry
  ``renewable_carriers``. Before this was done by commenting/uncommenting
  sub-sections in the ``renewable`` config section.

* Now, all carriers that should be extendable have to be listed in the config
  entry ``extendable_carriers``. Before, renewable carriers were always set to
  be extendable. For backwards compatibility, the workflow is still looking at
  the listed carriers under the ``renewable`` key. In the future, all of them
  have to be listed under ``extendable_carriers``.

* It is now possible to set conventional power plants as extendable by adding
  them to the list of extendable ``Generator`` carriers in the config.

* Listing conventional carriers in ``extendable_carriers`` but not in
  ``conventional_carriers``, sets the corresponding conventional power plants as
  extendable without a lower capacity bound of today's capacities.

* Now, conventional carriers have an assigned capital cost by default.

* The ``build_year`` and ``lifetime`` column are now defined for conventional
  power plants.

* Use updated SARAH-2 and ERA5 cutouts with slightly wider scope to east and
  additional variables.

* Resource definitions for memory usage now follow `Snakemake standard resource
  definition
  <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#standard-resources>`_
  ``mem_mb`` rather than ``mem``.

* The powerplants that have been shut down by 2021 are filtered out.

* Updated historical `EIA hydro generation data <https://www.eia.gov/international/data/world>`_.

* Network building is made deterministic by supplying a fixed random state to
  network clustering routines.

* Clustering strategies for generator and bus attributes can now be specified directly in the ``config/config.yaml``.

* Iterative solving with impedance updates is skipped if there are no expandable
  lines.

* The unused argument ``simple_hvdc_costs`` in :mod:`add_electricity` was
  removed.

* Switch from Germany to Belgium for continuous integration and tutorial to save
  resources.

* It is now possible to skip the progressbar for land eligibility calculations for additional speedup.

**Bugs and Compatibility**

* Fix crs bug. Change crs 4236 to 4326.

* ``powerplantmatching>=0.5.1`` is now required for ``IRENASTATS``.

* Update rasterio version to correctly calculate exclusion raster.

* It is now possible to run the workflow with only landlocked countries.

* Bugfixes for manual load adjustments across years.

* Enable parallel computing with new dask version.

* Restore compatibility of ``mock_snakemake`` with latest Snakemake versions.

* Script ``build_bus_regions``: move voronoi partition from vresutils to script.

* Script ``add_electricity``: remove ``vresutils.costdata.annuity`` dependency.

* Fix the plot_network snakemake rule.

* Compatibility with pandas 1.4. Address deprecations.

* Restore Windows compatibility by using ``shutil.move`` rather than ``mv``.


Synchronisation Release - Ukraine and Moldova (17th March 2022)
---------------------------------------------------------------

On March 16, 2022, the transmission networks of Ukraine and Moldova have
successfully been `synchronised with the continental European grid <https://www.entsoe.eu/news/2022/03/16/continental-europe-successful-synchronisation-with-ukraine-and-moldova-power-systems/>`_. We have taken
this as an opportunity to add the power systems of Ukraine and Moldova to
PyPSA-Eur. This includes:

.. image:: img/synchronisation.png
  :width: 500

* the transmission network topology from the `ENTSO-E interactive map <https://www.entsoe.eu/data/map/>`_.

* existing power plants (incl. nuclear, coal, gas and hydro) from the `powerplantmatching <https://github.com/fresna/powerplantmatching>`_ tool

* country-level load time series from ENTSO-E through the `OPSD platform <https://data.open-power-system-data.org/time_series/2020-10-06>`_, which are then distributed heuristically to substations by GDP and population density.

* wind and solar profiles based on ERA5 and SARAH-2 weather data

* hydro profiles based on historical `EIA generation data <https://www.eia.gov/international/data/world>`_

* a simplified calculation of wind and solar potentials based on the `Copernicus Land Cover dataset <https://land.copernicus.eu/global/products/lc>`_.

* electrical characteristics of 750 kV transmission lines

The Crimean power system is currently disconnected from the main Ukrainian grid and, hence, not included.

This release is not on the ``master`` branch. It can be used with

.. code-block:: bash

  git clone https://github.com/pypsa/pypsa-eur
  git checkout synchronisation-release


PyPSA-Eur 0.4.0 (22th September 2021)
-------------------------------------

**New Features and Changes**

* With this release, we change the license from copyleft GPLv3 to the more
  liberal MIT license with the consent of all contributors
  [`#276 <https://github.com/PyPSA/pypsa-eur/pull/276>`_].

* Switch to the new major ``atlite`` release v0.2.  The version upgrade comes
  along with significant speed up for the rule ``build_renewable_profiles.py``
  (~factor 2). A lot of the code which calculated the land-use availability is now
  outsourced and does not rely on ``glaes``, ``geokit`` anymore. This facilitates
  the environment building and version compatibility of ``gdal``, ``libgdal`` with
  other packages [`#224 <https://github.com/PyPSA/pypsa-eur/pull/224>`_].

* Implemented changes to ``n.snapshot_weightings`` in new PyPSA version v0.18
  (cf. `PyPSA/PyPSA/#227 <https://github.com/PyPSA/PyPSA/pull/227>`_)
  [`#259 <https://github.com/PyPSA/pypsa-eur/pull/259>`_].

* Add option to pre-aggregate nodes without power injections (positive or
  negative, i.e. generation or demand) to electrically closest nodes or neighbors
  in ``simplify_network``. Defaults to ``False``. This affects nodes that are no
  substations or have no offshore connection.

* In :mod:`simplify_network`, bus columns with no longer correct entries are
  removed (symbol, tags, under_construction, substation_lv, substation_off)
  [`#219 <https://github.com/PyPSA/pypsa-eur/pull/219>`_]

* Add option to include marginal costs of links representing fuel cells,
  electrolysis, and battery inverters
  [`#232 <https://github.com/PyPSA/pypsa-eur/pull/232>`_].

* The rule and script ``build_country_flh`` are removed as they are no longer
  used or maintained.

* The connection cost of generators in :mod:`simplify_network` are now reported
  in ``resources/connection_costs_s{simpl}.csv``
  [`#261 <https://github.com/PyPSA/pypsa-eur/pull/261>`_].

* The tutorial cutout was renamed from ``cutouts/europe-2013-era5.nc`` to
  ``cutouts/be-03-2013-era5.nc`` to accommodate tutorial and productive
  cutouts side-by-side.

* The flag ``keep_all_available_areas`` in the configuration for renewable
  potentials was deprecated and now defaults to ``True``.

* Update dependencies in ``envs/environment.yaml``
  [`#257 <https://github.com/PyPSA/pypsa-eur/pull/257>`_]

* Continuous integration testing switches to Github Actions from Travis CI
  [`#252 <https://github.com/PyPSA/pypsa-eur/pull/252>`_].

* Documentation on readthedocs.io is now built with ``pip`` only and no longer
  requires ``conda`` [`#267 <https://github.com/PyPSA/pypsa-eur/pull/267>`_].

* Use ``Citation.cff`` [`#273 <https://github.com/PyPSA/pypsa-eur/pull/273>`_].

**Bugs and Compatibility**


* Support for PyPSA v0.18 [`#268 <https://github.com/PyPSA/pypsa-eur/pull/268>`_].

* Minimum Python version set to ``3.8``.

* Removed ``six`` dependency [`#245 <https://github.com/PyPSA/pypsa-eur/pull/245>`_].

* Update :mod:`plot_network` and :mod:`make_summary` rules to latest PyPSA
  versions  [`#270 <https://github.com/PyPSA/pypsa-eur/pull/270>`_].

* Keep converter links to store components when using the ``ATK``
  wildcard and only remove DC links [`#214 <https://github.com/PyPSA/pypsa-eur/pull/214>`_].

* Value for ``co2base`` in ``config.yaml`` adjusted to 1.487e9 t CO2-eq
  (from 3.1e9 t CO2-eq). The new value represents emissions related to the
  electricity sector for EU+UK+Balkan. The old value was too high and used when
  the emissions wildcard in ``{opts}`` was used
  [`#233 <https://github.com/PyPSA/pypsa-eur/pull/233>`_].

* Add escape in :mod:`base_network` if all TYNDP links are already
  contained in the network
  [`#246 <https://github.com/PyPSA/pypsa-eur/pull/246>`_].

* In :mod:`solve_operations_network` the optimised capacities are now
  fixed for all extendable links, not only HVDC links
  [`#244 <https://github.com/PyPSA/pypsa-eur/pull/244>`_].

* The ``focus_weights`` are now also considered when pre-clustering in
  the :mod:`simplify_network` rule
  [`#241 <https://github.com/PyPSA/pypsa-eur/pull/241>`_].

* in :mod:`build_renewable_profile` where offshore wind profiles could
  no longer be created [`#249 <https://github.com/PyPSA/pypsa-eur/pull/249>`_].

* Lower expansion limit of extendable carriers is now set to the
  existing capacity, i.e. ``p_nom_min = p_nom`` (0 before). Simultaneously, the
  upper limit (``p_nom_max``) is now the maximum of the installed capacity
  (``p_nom``) and the previous estimate based on land availability (``p_nom_max``)
  [`#260 <https://github.com/PyPSA/pypsa-eur/pull/260>`_].

* Solving an operations network now includes optimized store capacities
  as well. Before only lines, links, generators and storage units were considered
  [`#269 <https://github.com/PyPSA/pypsa-eur/pull/269>`_].

* With ``load_shedding: true`` in the solving options of ``config.yaml``
  load shedding generators are only added at the AC buses, excluding buses for H2
  and battery stores [`#269 <https://github.com/PyPSA/pypsa-eur/pull/269>`_].

* Delete duplicated capital costs at battery discharge link
  [`#240 <https://github.com/PyPSA/pypsa-eur/pull/240>`_].

* Propagate the solver log file name to the solver. Previously, the
  PyPSA network solving functions were not told about the solver logfile specified
  in the Snakemake file [`#247 <https://github.com/PyPSA/pypsa-eur/pull/247>`_]

PyPSA-Eur 0.3.0 (7th December 2020)
-----------------------------------

**New Features**

Using the ``{opts}`` wildcard for scenario:

* An option is introduced which adds constraints such that each country or node produces on average a minimal share of its total consumption itself.
  For example ``EQ0.5c`` set in the ``{opts}`` wildcard requires each country to produce on average at least 50% of its consumption. Additionally,
  the option ``ATK`` requires autarky at each node and removes all means of power transmission through lines and links. ``ATKc`` only removes
  cross-border transfer capacities.
  [`#166 <https://github.com/PyPSA/pypsa-eur/pull/166>`_].

* Added an option to alter the capital cost (``c``) or installable potentials (``p``) of carriers by a factor via ``carrier+{c,p}factor`` in the ``{opts}`` wildcard.
  This can be useful for exploring uncertain cost parameters.
  Example: ``solar+c0.5`` reduces the capital cost of solar to 50% of original values
  [`#167 <https://github.com/PyPSA/pypsa-eur/pull/167>`_, `#207 <https://github.com/PyPSA/pypsa-eur/pull/207>`_].

* Added an option to the ``{opts}`` wildcard that applies a time series segmentation algorithm based on renewables, hydro inflow and load time series
  to produce a given total number of adjacent snapshots of varying lengths.
  This feature is an alternative to downsampling the temporal resolution by simply averaging and
  uses the `tsam <https://tsam.readthedocs.io/en/latest/index.html>`_ package
  [`#186 <https://github.com/PyPSA/pypsa-eur/pull/186>`_].


More OPSD integration:

* Add renewable power plants from `OPSD <https://data.open-power-system-data.org/renewable_power_plants/2020-08-25>`_ to the network for specified technologies.
  This will overwrite the capacities calculated from the heuristic approach in :func:`estimate_renewable_capacities()`
  [`#212 <https://github.com/PyPSA/pypsa-eur/pull/212>`_].

* Electricity consumption data is now retrieved directly from the `OPSD website <https://data.open-power-system-data.org/time_series/2019-06-05>`_ using the rule :mod:`build_electricity_demand`.
  The user can decide whether to take the ENTSO-E power statistics data (default) or the ENTSO-E transparency data
  [`#211 <https://github.com/PyPSA/pypsa-eur/pull/211>`_].

Other:

* Added an option to use custom busmaps in rule :mod:`cluster_network`. To use this feature set ``enable: custom_busmap: true``.
  Then, the rule looks for custom busmaps at ``data/custom_busmap_elec_s{simpl}_{clusters}.csv``,
  which should have the same format as ``resources/busmap_elec_s{simpl}_{clusters}.csv``.
  i.e. the index should contain the buses of ``networks/elec_s{simpl}.nc``
  [`#193 <https://github.com/PyPSA/pypsa-eur/pull/193>`_].

* Line and link capacities can be capped in the ``config.yaml`` at ``lines: s_nom_max:`` and ``links: p_nom_max``:
  [`#166 <https://github.com/PyPSA/pypsa-eur/pull/166>`_].

* Added Google Cloud Platform tutorial (for Windows users)
  [`#177 <https://github.com/PyPSA/pypsa-eur/pull/177>`_].

**Changes**

* Don't remove capital costs from lines and links, when imposing a line volume limit (``lv``) or a line cost limit (``lc``).
  Previously, these were removed to move the expansion in direction of the limit
  [`#183 <https://github.com/PyPSA/pypsa-eur/pull/183>`_].

* The mappings for clustered lines and buses produced by the :mod:`simplify_network` and :mod:`cluster_network` rules
  changed from Hierarchical Data Format (``.h5``) to Comma-Separated Values format (``.csv``) for ease of use.
  [`#198 <https://github.com/PyPSA/pypsa-eur/pull/198>`_]

* The N-1 security margin for transmission lines is now fixed to a provided value in ``config.yaml``,
  removing an undocumented linear interpolation between 0.5 and 0.7 in the range between 37 and 200 nodes.
  [`#199 <https://github.com/PyPSA/pypsa-eur/pull/199>`_].

* Modelling hydrogen and battery storage with Store and Link components is now the default,
  rather than using StorageUnit components with fixed power-to-energy ratio
  [`#205 <https://github.com/PyPSA/pypsa-eur/pull/205>`_].

* Use ``mamba`` (https://github.com/mamba-org/mamba) for faster Travis CI builds
  [`#196 <https://github.com/PyPSA/pypsa-eur/pull/196>`_].

* Multiple smaller changes: Removed unused ``{network}`` wildcard, moved environment files to dedicated ``envs`` folder,
  removed sector-coupling components from configuration files, updated documentation colors, minor refactoring and code cleaning
  [`#190 <https://github.com/PyPSA/pypsa-eur/pull 190>`_].

**Bugs and Compatibility**

* Add compatibility for pyomo 5.7.0 in :mod:`cluster_network` and :mod:`simplify_network`
  [`#172 <https://github.com/PyPSA/pypsa-eur/pull/172>`_].

* Fixed a bug for storage units such that individual store and dispatch efficiencies are correctly taken account of rather than only their round-trip efficiencies.
  In the cost database (``data/costs.csv``) the efficiency of battery inverters should be stated as per discharge/charge rather than per roundtrip
  [`#202 <https://github.com/PyPSA/pypsa-eur/pull/202>`_].

* Corrected exogenous emission price setting (in ``config: cost: emission price:``),
  which now correctly accounts for the efficiency and effective emission of the generators
  [`#171 <https://github.com/PyPSA/pypsa-eur/pull/171>`_].

* Corrected HVDC link connections (a) between Norway and Denmark and (b) mainland Italy, Corsica (FR) and Sardinia (IT)
  as well as for East-Western and Anglo-Scottish interconnectors
  [`#181 <https://github.com/PyPSA/pypsa-eur/pull/181>`_, `#206 <https://github.com/PyPSA/pypsa-eur/pull/206>`_].

* Fix bug of clustering ``offwind-{ac,dc}`` generators in the option of high-resolution generators for renewables.
  Now, there are more sites for ``offwind-{ac,dc}`` available than network nodes.
  Before, they were clustered to the resolution of the network (``elec_s1024_37m.nc``: 37 network nodes, 1024 generators)
  [`#191 <https://github.com/PyPSA/pypsa-eur/pull/191>`_].

* Raise a warning if ``tech_colors`` in the config are not defined for all carriers
  [`#178 <https://github.com/PyPSA/pypsa-eur/pull/178>`_].


PyPSA-Eur 0.2.0 (8th June 2020)
-------------------------------

* The optimization is now performed using the ``pyomo=False`` setting in the :func:`pypsa.lopf.network_lopf`. This speeds up the solving process significantly and consumes much less memory. The inclusion of additional constraints were adjusted to the new implementation. They are all passed to the :func:`network_lopf` function via the ``extra_functionality`` argument. The rule ``trace_solve_network`` was integrated into the rule :mod:`solve_network` and can be activated via configuration with ``solving: options: track_iterations: true``. The charging and discharging capacities of batteries modelled as store-link combination are now coupled [`#116 <https://github.com/PyPSA/pypsa-eur/pull/116>`_].

* An updated extract of the `ENTSO-E Transmission System Map <https://www.entsoe.eu/data/map/>`_ (including Malta) was added to the repository using the `GridKit <https://github.com/PyPSA/GridKit>`_ tool. This tool has been updated to retrieve up-to-date map extracts using a single `script <https://github.com/PyPSA/GridKit/blob/master/entsoe/runall_in_docker.sh>`_. The update extract features 5322 buses, 6574 lines, 46 links. [`#118 <https://github.com/PyPSA/pypsa-eur/pull/118>`_].

* Added `FSFE REUSE <https://reuse.software>`_ compliant license information. Documentation now licensed under CC-BY-4.0 [`#160 <https://github.com/PyPSA/pypsa-eur/pull/160>`_].

* Added a 30 minute `video introduction <https://pypsa-eur.readthedocs.io/en/latest/introduction.html>`_ and a 20 minute `video tutorial <https://pypsa-eur.readthedocs.io/en/latest/tutorial.html>`_

* Networks now store a color and a nicely formatted name for each carrier, accessible via ``n.carrier['color']`` and ``n.carrier['nice_name'] ``(networks after ``elec.nc``).

* Added an option to skip iterative solving usually performed to update the line impedances of expanded lines at ``solving: options: skip_iterations:``.

* ``snakemake`` rules for retrieving cutouts and the natura raster can now be disabled independently from their respective rules to build them; via ``config.*yaml`` [`#136 <https://github.com/PyPSA/pypsa-eur/pull/136>`_].

* Removed the ``id`` column for custom power plants in ``data/custom_powerplants.csv`` to avoid custom power plants with conflicting ids getting attached to the wrong bus [`#131 <https://github.com/PyPSA/pypsa-eur/pull/131>`_].

* Add option ``renewables: {carrier}: keep_all_available_areas:`` to use all available weather cells for renewable profile and potential generation. The default ignores weather cells where only less than 1 MW can be installed  [`#150 <https://github.com/PyPSA/pypsa-eur/pull/150>`_].

* Added a function ``_helpers.load_network()`` which loads a network with overridden components specified in ``snakemake.config['override_components']`` [`#128 <https://github.com/PyPSA/pypsa-eur/pull/128>`_].

* Bugfix in  :mod:`base_network` which now finds all closest links, not only the first entry [`#143 <https://github.com/PyPSA/pypsa-eur/pull/143>`_].

* Bugfix in :mod:`cluster_network` which now skips recalculation of link parameters if there are no links  [`#149 <https://github.com/PyPSA/pypsa-eur/pull/149>`_].

* Added information on pull requests to contribution guidelines [`#151 <https://github.com/PyPSA/pypsa-eur/pull/151>`_].

* Improved documentation on open-source solver setup and added usage warnings.

* Updated ``conda`` environment regarding ``pypsa``, ``pyproj``, ``gurobi``, ``lxml``. This release requires PyPSA v0.17.0.

PyPSA-Eur 0.1.0 (9th January 2020)
----------------------------------

This is the first release of PyPSA-Eur, a model of the European power system at the transmission network level. Recent changes include:

* Documentation on installation, workflows and configuration settings is now available online at `pypsa-eur.readthedocs.io <pypsa-eur.readthedocs.io>`_ [`#65 <https://github.com/PyPSA/pypsa-eur/pull/65>`_].

* The ``conda`` environment files were updated and extended [`#81 <https://github.com/PyPSA/pypsa-eur/pull/81>`_].

* The power plant database was updated with extensive filtering options via ``pandas.query`` functionality [`#84 <https://github.com/PyPSA/pypsa-eur/pull/84>`_ and `#94 <https://github.com/PyPSA/pypsa-eur/pull/94>`_].

* Continuous integration testing with `Travis CI <https://travis-ci.org>`_ is now included for Linux, Mac and Windows [`#82 <https://github.com/PyPSA/pypsa-eur/pull/82>`_].

* Data dependencies were moved to `zenodo <https://zenodo.org/>`_ and are now versioned [`#60 <https://github.com/PyPSA/pypsa-eur/issues/60>`_].

* Data dependencies are now retrieved directly from within the snakemake workflow [`#86 <https://github.com/PyPSA/pypsa-eur/pull/86>`_].

* Emission prices can be added to marginal costs of generators through the keywords ``Ep`` in the ``{opts}`` wildcard [`#100 <https://github.com/PyPSA/pypsa-eur/pull/100>`_].

* An option is introduced to add extendable nuclear power plants to the network [`#98 <https://github.com/PyPSA/pypsa-eur/pull/98>`_].

* Focus weights can now be specified for particular countries for the network clustering, which allows to set a proportion of the total number of clusters for particular countries [`#87 <https://github.com/PyPSA/pypsa-eur/pull/87>`_].

* A new rule :mod:`add_extra_components` allows to add additional components to the network only after clustering. It is thereby possible to model storage units (e.g. battery and hydrogen) in more detail via a combination of ``Store``, ``Link`` and ``Bus`` elements [`#97 <https://github.com/PyPSA/pypsa-eur/pull/97>`_].

* Hydrogen pipelines (including cost assumptions) can now be added alongside clustered network connections in the rule :mod:`add_extra_components` . Set ``electricity: extendable_carriers: Link: [H2 pipeline]`` and ensure hydrogen storage is modelled as a ``Store``. This is a first simplified stage [`#108 <https://github.com/PyPSA/pypsa-eur/pull/108>`_].

* Logfiles for all rules of the ``snakemake`` workflow are now written in the folder ``log/`` [`#102 <https://github.com/PyPSA/pypsa-eur/pull/102>`_].

* The new function ``_helpers.mock_snakemake`` creates a ``snakemake`` object which mimics the actual ``snakemake`` object produced by workflow by parsing the ``Snakefile`` and setting all paths for inputs, outputs, and logs. This allows running all scripts within a (I)python terminal (or just by calling ``python <script-name>``) and thereby facilitates developing and debugging scripts significantly [`#107 <https://github.com/PyPSA/pypsa-eur/pull/107>`_].


PyPSA-Eur-Sec Releases (pre-merge)
==================================

PyPSA-Eur-Sec 0.7.0 (16th February 2023)
----------------------------------------

This release includes many new features. Highlights include new gas
infrastructure data with retrofitting options for hydrogen transport, improved
carbon management and infrastructure planning, regionalised potentials for
hydrogen underground storage and carbon sequestration, new applications for
biomass, and explicit modelling of methanol and ammonia as separate energy
carriers.

This release is known to work with `PyPSA-Eur
<https://github.com/PyPSA/pypsa-eur>`_ Version 0.7.0 and `Technology Data
<https://github.com/PyPSA/technology-data>`_ Version 0.5.0.

**Gas Transmission Network**

* New rule ``retrieve_gas_infrastructure_data`` that downloads and extracts the
  SciGRID_gas `IGGIELGN <https://zenodo.org/record/4767098>`_ dataset from
  zenodo. It includes data on the transmission routes, pipe diameters,
  capacities, pressure, and whether the pipeline is bidirectional and carries
  H-Gas or L-Gas.

* New rule ``build_gas_network`` processes and cleans the pipeline data from
  SciGRID_gas. Missing or uncertain pipeline capacities can be inferred by
  diameter.

* New rule ``build_gas_input_locations`` compiles the LNG import capacities
  (from the Global Energy Monitor's `Europe Gas Tracker
  <https://globalenergymonitor.org/projects/europe-gas-tracker/>`_, pipeline
  entry capacities and local production capacities for each region of the model.
  These are the regions where fossil gas can eventually enter the model.

* New rule ``cluster_gas_network`` that clusters the gas transmission network
  data to the model resolution. Cross-regional pipeline capacities are
  aggregated (while pressure and diameter compatibility is ignored),
  intra-regional pipelines are dropped. Lengths are recalculated based on the
  regions' centroids.

* With the option ``sector: gas_network:``, the existing gas network is added
  with a lossless transport model. A length-weighted `k-edge augmentation
  algorithm
  <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation>`_
  can be run to add new candidate gas pipelines such that all regions of the
  model can be connected to the gas network. The number of candidates can be
  controlled via the setting ``sector: gas_network_connectivity_upgrade:``. When
  the gas network is activated, all the gas demands are regionally disaggregated
  as well.

* New constraint allows endogenous retrofitting of gas pipelines to hydrogen
  pipelines. This option is activated via the setting ``sector: H2_retrofit:``.
  For every unit of gas pipeline capacity dismantled, ``sector:
  H2_retrofit_capacity_per_CH4`` units are made available as hydrogen pipeline
  capacity in the corresponding corridor. These repurposed hydrogen pipelines
  have lower costs than new hydrogen pipelines. Both new and repurposed
  pipelines can be built simultaneously. The retrofitting option ``sector:
  H2_retrofit:`` also works with a copperplated methane infrastructure, i.e.
  when ``sector: gas_network: false``.

* New hydrogen pipelines can now be built where there are already power or gas
  transmission routes. Previously, only the electricity transmission routes were
  considered.

**Carbon Management and Biomass**

* Add option to spatially resolve carrier representing stored carbon dioxide
  (``co2_spatial``). This allows for more detailed modelling of CCUTS, e.g.
  regarding the capturing of industrial process emissions, usage as feedstock
  for electrofuels, transport of carbon dioxide, and geological sequestration
  sites.

* Add option for regionally-resolved geological carbon dioxide sequestration
  potentials through new rule ``build_sequestration_potentials`` based on
  `CO2StoP <https://setis.ec.europa.eu/european-co2-storage-database_en>`_. This
  can be controlled in the section ``regional_co2_sequestration_potential`` of
  the ``config.yaml``. It includes options to select the level of conservatism,
  whether onshore potentials should be included, the respective upper and lower
  limits per region, and an annualisation parameter for the cumulative
  potential. The defaults are preliminary and will be validated the next
  release.

* Add option to sweep the global CO2 sequestration potentials with keyword
  ``seq200`` in the ``{sector_opts}`` wildcard (for limit of 200 Mt CO2).

* Add option to include `Allam cycle gas power plants
  <https://en.wikipedia.org/wiki/Allam_power_cycle>`_ (``allam_cycle``).

* Add option for planning a new carbon dioxide network (``co2network``).

* Separate option to regionally resolve biomass (``biomass_spatial``) from
  option to allow biomass transport (``biomass_transport``).

* Add option for biomass boilers (wood pellets) for decentral heating.

* Add option for BioSNG (methane from biomass) with and without carbon capture.

* Add option for BtL (biomass to liquid fuel/oil) with and without carbon
  capture.


**Other new features**

* Add regionalised hydrogen salt cavern storage potentials from `Technical
  Potential of Salt Caverns for Hydrogen Storage in Europe
  <https://doi.org/10.20944/preprints201910.0187.v1>`_. This data is compiled in
  a new rule ``build_salt_cavern_potentials``.

* Add option to resolve ammonia as separate energy carrier with Haber-Bosch
  synthesis, ammonia cracking, storage and industrial demand. The ammonia
  carrier can be nodally resolved or copperplated across Europe (see
  ``ammonia``).

* Add methanol as energy carrier, methanolisation as process, and option for
  methanol demand in shipping sector.

* Shipping demand now defaults to methanol rather than liquefied hydrogen
  until 2050.

* Demand for liquid hydrogen in international shipping is now geographically
  distributed by port trade volumes in a new rule ``build_shipping_demand``
  using data from the `World Bank Data Catalogue
  <https://datacatalog.worldbank.org/search/dataset/0038118/Global---International-Ports>`_.
  Domestic shipping remains distributed by population.

* Add option to aggregate network temporally using representative snapshots or
  segments (with `tsam <https://github.com/FZJ-IEK3-VSA/tsam>`_).

* Add option for minimum part load for Fischer-Tropsch plants (default: 90%) and
  methanolisation plants (default: 50%).

* Add option to use waste heat of electrolysis in district heating networks
  (``use_electrolysis_waste_heat``).

* Add option for coal CHPs with carbon capture (see ``coal_cc``).

* In overnight optimisation, it is now possible to specify a year for the
  technology cost projections separate from the planning horizon.

* New config options for changing energy demands in aviation
  (``aviation_demand_factor``) and HVC industry (``HVC_demand_factor``), as well
  as explicit ICE shares for land transport (``land_transport_ice_share``) and
  agriculture machinery (``agriculture_machinery_oil_share``).

* It is now possible to merge residential and services heat buses to reduce the
  problem size (see ``cluster_heat_nodes``).

* Added option to tweak (almost) any configuration parameter through the
  ``{sector_opts}`` wildcard. The regional_co2_sequestration_potential is
  triggered by the prefix ``CF+`` after which it is possible to pipe to any
  setting that does not contain underscores (``_``). Example:
  ``CF+sector+v2g+false`` disables vehicle-to-grid flexibility.

* Option ``retrieve_sector_databundle`` to automatically retrieve and extract
  data bundle.

* Removed the need to clone ``technology-data`` repository in a parallel
  directory. The new approach automatically retrieves the technology data from
  remote in the rule ``retrieve_cost_data``.

* Improved network plots including better legends, hydrogen retrofitting network
  display, and change to EqualEarth projection. A new color scheme for
  technologies was also introduced.

* Add two new rules ``build_transport_demand`` and
  ``build_population_weighted_energy_totals`` using code previously contained in
  ``prepare_sector_network``.

* Rules that convert weather data with ``atlite`` now largely run separately for
  categories residential, rural and total.

* Units are assigned to the buses. These only provide a better understanding.
  The specifications of the units are not taken into account in the
  optimisation, which means that no automatic conversion of units takes place.

* Configuration file and wildcards are now stored under ``n.meta`` in every
  PyPSA network.

* Updated `data bundle
  <https://zenodo.org/record/5824485/files/pypsa-eur-sec-data-bundle.tar.gz>`_
  that includes the hydrogan salt cavern storage potentials.

* Updated and extended documentation in
  <https://pypsa-eur-sec.readthedocs.io/en/latest/>

* Added new rule ``copy_conda_env`` that exports a list of packages with which
  the workflow was executed.

* Add basic continuous integration using Github Actions.

* Add basic ``rsync`` setup.

**Bugfixes**

* The CO2 sequestration limit implemented as GlobalConstraint (introduced in the
  previous version) caused a failure to read in the shadow prices of other
  global constraints.

* Correct capital cost of Fischer-Tropsch according to new units in
  ``technology-data`` repository.

* Fix unit conversion error for thermal energy storage.

* For myopic pathway optimisation, set optimised capacities of power grid
  expansion of previous iteration as minimum capacity for next iteration.

* Further rather minor bugfixes for myopic optimisation code (see `#256
  <https://github.com/PyPSA/pypsa-eur-sec/pull/256>`_).


Many thanks to all who contributed to this release!


PyPSA-Eur-Sec 0.6.0 (4 October 2021)
------------------------------------

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
-----------------------------------

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
----------------------------------------

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
-----------------------------------------

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
--------------------------------------

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
-----------------------------------

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
----------------------------------------

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
---------------------------------------

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

* Checkout a new release branch ``git checkout -b release-v0.x.x``.

* Finalise release notes at ``doc/release_notes.rst``.

* Update ``envs/environment.fixed.yaml`` via
  ``conda env export -n pypsa-eur -f envs/environment.fixed.yaml --no-builds``
  from an up-to-date ``pypsa-eur`` environment.

* Update version number in ``doc/conf.py``, ``CITATION.cff`` and ``*config.*.yaml``.

* Make a ``git commit``.

* Open, review and merge pull request for branch ``release-v0.x.x``.
  Make sure to close issues and PRs or the release milestone with it (e.g. closes #X).

* Tag a release on Github via ``git tag v0.x.x``, ``git push``, ``git push --tags``. Include release notes in the tag message.

* Make a `GitHub release <https://github.com/PyPSA/pypsa-eur-sec/releases>`_, which automatically triggers archiving to the `zenodo code repository <https://doi.org/10.5281/zenodo.3520874>`_ with `MIT license <https://opensource.org/licenses/MIT>`_.

* Send announcement on the `PyPSA mailing list <https://groups.google.com/forum/#!forum/pypsa>`_.
