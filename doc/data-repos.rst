..
  SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

###########
Repository
###########

Data in this section is included in the PyPSA-Eur repository in the ``data`` folder.

``data/entsoegridkit``

- **Source:** ENTSO-E
- **Link:** https://www.entsoe.eu/data/map/, extracted with https://github.com/PyPSA/GridKit/tree/master/entsoe
- **License:** unknown
- **Description:** Data of high-voltage transmission grid in Europe from ENTSO-E.

``data/existing_infrastructure``

- **Source:** European Commission DG ENER; Mapping and analyses of the current and future (2020 - 2030) heating/cooling fuel deployment
- **Link:** https://energy.ec.europa.eu/publications/mapping-and-analyses-current-and-future-2020-2030-heatingcooling-fuel-deployment-fossilrenewables-1_en
- **License:** CC-BY 4.0 (`reference <https://commission.europa.eu/legal-notice_en>`__)
- **Description:** Contains country-level data on existing heating infrastructure, i.e. gas, oil, coal boilers, resistive heaters, air- and ground-sourced heat pumps.

``data/retro/comparative_level_investment.csv``

- **Source:** Eurostat
- **Link:** https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Comparative_price_levels_for_investment
- **License:** `custom <https://ec.europa.eu/eurostat/about-us/policies/copyright>`__
- **Description:** Contains data on comparative price levels for investment in Europe.

``data/retro/data_building_stock.csv``

- **Source:** Simon Pezzutto, Stefano Zambotti, Silvia Croce, Pietro Zambelli,
  Giulia Garegnani, Chiara Scaramuzzino, Ramón Pascual Pascuas, Alyona
  Zubaryeva, Franziska Haas, Dagmar Exner (EURAC), Andreas Müller (e-think),
  Michael Hartner (TUW), Tobias Fleiter, Anna-Lena Klingler, Matthias Kühnbach,
  Pia Manz, Simon Marwitz, Matthias Rehfeldt, Jan Steinbach, Eftim Popovski
  (Fraunhofer ISI) Reviewed by Lukas Kranzl, Sara Fritz (TUW) Hotmaps Project,
  D2.3 WP2 Report - Open Data Set for the EU28, 2018 www.hotmaps-project.eu
- **Link:** https://gitlab.com/hotmaps/building-stock
- **License:** CC-BY 4.0
- **Description:** Contains data on European building stock.

``data/retro/electricity_taxes_eu.csv``

- **Source:** Eurostat
- **Link:** https://ec.europa.eu/eurostat/databrowser/view/NRG_PC_204/default/table?lang=en
- **License:** `custom <https://ec.europa.eu/eurostat/about-us/policies/copyright>`__
- **Description:** Electricity prices for household consumers.

``data/retro/{floor_area_missing,u_values_poland}.csv``

- **Source:** EU Building Stock Observatory
- **Link:** https://data.europa.eu/euodp/de/data/dataset/building-stock-observatory
- **License:** `custom <https://data.europa.eu/data/datasets/building-stock-observatory?locale=en>`__
- **Description:** The EU Building Stock Observatory monitors the energy
  performance of buildings across Europe. It assesses improvements in the energy
  efficiency of buildings and the impact of this on the actual energy
  consumption of the buildings sector overall.

``data/retro/retro_cost_germany.csv``

- **Source:** Institut Wohnen und Umwelt (IWU)
- **Link:** https://www.iwu.de/forschung/handlungslogiken/kosten-energierelevanter-bau-und-anlagenteile-bei-modernisierung/
- **License:** unknown
- **Description:** Contains thermal envelop costs for retrofitting buildings in
  Germany.

``data/retro/window_assumptions.csv``

- **Source:** ifeu, Fraunhofer IEE and Consentec (2018): Building sector
  Efficiency: A crucial Component of the Energy Transition. A study commissioned
  by Agora Energiewende.
- **Link:** https://www.agora-energiewende.de/en/publications/building-sector-efficiency-a-crucial-component-of-the-energy-transition/
- **License:** unknown
- **Description:** Contains data on physical parameters of double- and triple-glazed windows.

``data/transmission_projects/nep``

- **Source:** German Federal Network Agency (Bundesnetzagentur, BNetzA)
- **Link:** https://data.netzausbau.de/2037-2023/NEP/NEP_2037_2045_Bestaetigung.pdf
- **License:** unknown
- **Description:** Contains transmission projects in Europe from German network development plan (Netzentwicklungsplan).

``data/transmission_projects/tyndp2020``

- **Source:** ENTSO-E
- **Link:** https://tyndp2020-project-platform.azurewebsites.net/projectsheets
- **License:** unknown
- **Description:** Contains transmission projects in Europe from ENTSO-E Ten Year Network Development Plan (TYNDP).

``data/ammonia_plants.csv``

- **Source:** manually collected, mostly from ICIS
- **Link:** https://www.icis.com/explore/resources/news/2023/01/18/10846094/insight-poor-demand-high-costs-stifle-europe-industry-despite-falling-gas-prices/
- **License:** CC-BY 4.0 (for compiled dataset)
- **Description:** Locations and production capacities of ammonia plants in Europe.

``data/attributed_ports.json``

- **Source:** World Bank
- **Link:** https://datacatalog.worldbank.org/search/dataset/0038118/Global---International-Ports
- **License:** CC-BY 4.0 (`reference <https://datacatalog.worldbank.org/search/dataset/0038118/Global---International-Ports>`__)
- **Description:** International ports with attributes describing name, port functions, total capacity and location.

``data/cement_plants-noneu.csv``

- **Source:** manually collected, mostly from USGS
- **Link:** https://www.usgs.gov/centers/national-minerals-information-center/international-minerals-statistics-and-information
- **License:** CC0 (`reference <https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits>`__)
- **Description:** Contains energy balances for Europe.

``data/ch_cantons.csv``

- **Source:** Wikipedia
- **Link:** https://en.wikipedia.org/wiki/Data_codes_for_Switzerland
- **License:** CC-BY-SA 4.0
- **Description:** Contains NUTS codes for regions in Switzerland.

``data/ch_industrial_production_per_subsector.csv``

- **Source:** Swiss Federal Office of Energy (SFOE)
- **Link:** https://pubdb.bfe.admin.ch/de/publication/download/11817
- **License:** `custom <https://www.admin.ch/gov/de/start/rechtliches.html>`__
- **Description:** Contains energy consumption in industry and the service sector in Switzerland.

``data/district_heat_share.csv``

- **Source:** Euroheat & Power
- **Link:** https://www.euroheat.org/knowledge-hub/country-profiles
- **License:** unknown
- **Description:** Contains district heating shares for European countries.

``data/egs_costs.json``

- **Source:** Arman Aghahosseini, Christian Breyer, From hot rock to useful
  energy: A global estimate of enhanced geothermal systems potential, Applied
  Energy, Volume 279, 2020, 115769.
- **Link:** https://doi.org/10.1016/j.apenergy.2020.115769
- **License:** unknown
- **Description:** Contains rastered potentials and capital costs for enhanced geothermal electricity generation in Europe.

``data/eia_hydro_annual_capacity.csv``

- **Source:** Energy Information Agency (EIA)
- **Link:** https://www.eia.gov/international/data/world/electricity/electricity-generation
- **License:** CC0 (`reference <https://www.eia.gov/about/copyrights_reuse.php>`__)
- **Description:** Contains country-level hydro-electric capacity for Europe by year.

``data/eia_hydro_annual_generation.csv``

- **Source:** Energy Information Agency (EIA)
- **Link:** https://www.eia.gov/international/data/world/electricity/electricity-generation
- **License:** CC0 (`reference <https://www.eia.gov/about/copyrights_reuse.php>`__)
- **Description:** Contains country-level hydro-electric generato for Europe by year.

``data/era5-annual-HDD-per-country.csv``

- **Source:** Neumann, Fabian
- **Link:** https://gist.github.com/fneum/d99e24e19da423038fd55fe3a4ddf875
- **License:** CC-BY 4.0
- **Description:** Contains country-level annual sum of heating degree days in
  Europe. Used for rescaling heat demand in weather years not covered by energy
  balance statistics.

``data/era5-annual-runoff-per-country.csv``

- **Source:** Neumann, Fabian
- **Link:** https://gist.github.com/fneum/d99e24e19da423038fd55fe3a4ddf875
- **License:** CC-BY 4.0
- **Description:** Contains country-level annual sum of runoff in Europe. Used
  for rescaling hydro-electricity availability in weather years not covered by
  EIA hydro-generation statistics.

``data/gr-e-11.03.02.01.01-cc.csv``

- **Source:** Swiss Federal Statistics Office
- **Link:** https://www.bfs.admin.ch/asset/de/30305426
- **License:** `custom (OPEN BY ASK) <https://www.bfs.admin.ch/bfs/en/home/fso/swiss-federal-statistical-office/terms-of-use.html>`__
- **Description:** Stock of road motor vehicles in Switzerland.

``data/heat_load_profile_BDEW.csv``

- **Source:** oemof/demandlib
- **Link:** https://github.com/oemof/demandlib
- **License:** unknown
- **Description:** Contains standard heat load profiles based on data from BDEW (German Association of Energy and Water Industries).

``data/hydro_capacities.csv``

.. warning::
   The provenance of the data is unclear. We will improve this in the future.

``data/links_p_nom.csv``

- **Source:** Wikipedia
- **Link:** https://en.wikipedia.org/wiki/List_of_HVDC_projects
- **License:** CC-BY-SA 4.0
- **Description:** Contains list of HVDC transmission line projects.

``data/nuclear_p_max_pu.csv``

- **Source:** International Atomic Energy Agency (IAEA)
- **Link:** https://pris.iaea.org/PRIS/WorldStatistics/ThreeYrsEnergyAvailabilityFactor.aspx
- **License:** `custom <https://www.iaea.org/about/terms-of-use>`__
- **Description:** Country-level nuclear power plant availability factors.

``data/refineries-noneu.csv``

- **Source:** manually collected, mostly from Energy Information Agency (EIA)
- **Link:** https://www.eia.gov/petroleum/refinerycapacity/table3.pdf
- **License:** CC0 (`reference <https://www.eia.gov/about/copyrights_reuse.php>`__)
- **Description:** Contains locations and capacities of oil refineries in Europe.

``data/switzerland-new_format-all_years.csv``

- **Source:** Swiss Federal Office of Energy (SFOE)
- **Link:** https://www.bfe.admin.ch/bfe/de/home/versorgung/statistik-und-geodaten/energiestatistiken/energieverbrauch-nach-verwendungszweck.html/
- **License:** `custom <https://www.admin.ch/gov/de/start/rechtliches.html>`__
- **Description:** Contains energy consumption by sector / application for Switzerland.

``data/unit_commitment.csv``

- **Source:** `DIW
  <https://www.diw.de/documents/publikationen/73/diw_01.c.424566.de/diw_datadoc_2013-068.pdf>`__,
  `Agora Energiewende
  <https://www.agora-energiewende.de/fileadmin/Projekte/2017/Flexibility_in_thermal_plants/115_flexibility-report-WEB.pdf>`__,
  `Schill et al. (2017)
  <https://static-content.springer.com/esm/art%3A10.1038%2Fnenergy.2017.50/MediaObjects/41560_2017_BFnenergy201750_MOESM196_ESM.pdf>`__,
  `Martin (2022) <https://zenodo.org/records/6421682>`__
- **Link:** https://github.com/lisazeyen/hourly_vs_annually/blob/b67ca9222711372d8ab6cd58f9ebe7bc637939bf/scripts/solve_network.py#L554
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/biomass_transport_costs_supply_chain{1,2}.csv``

- **Source:** European Commission Joint Research Centre (JRC)
- **Link:** https://publications.jrc.ec.europa.eu/repository/handle/JRC98626
- **License:** CC-BY 4.0 (`reference <https://commission.europa.eu/legal-notice_en#copyright-notice>`__)
- **Description:** Contains transport costs for different types of biomass.
