..
  SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Data Sources
##########################################

PyPSA-Eur is combiled from a variety of data sources. The following table provides an
overview of the data sources used in PyPSA-Eur. Different licenses apply to the
data sources.

Zenodo data bundle
=======================

Data in this section is downloaded and extracted from the Zenodo data bundle
(https://zenodo.org/records/12760663). Files included in the data bundle are too
large to be placed directly in the repository, have been reduced in spatial
scope to reduce file size, or are not provided through stable URLs elsewhere.

``data/bundle/je-e-21.03.02.xls``

- **Source:** Swiss Federal Statistics Office
- **Link:** https://www.bfs.admin.ch/bfs/en/home/news/whats-new.assetdetail.7786557.html#context-sidebar
- **License:**  `custom (OPEN BY ASK) <https://www.bfs.admin.ch/bfs/en/home/fso/swiss-federal-statistical-office/terms-of-use.html>`__
- **Description:** Population and GDP data for Swiss Cantons.

``data/bundle/NUTS_2013_60M_SH``

- **Source:** GISCO
- **Link:** https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/
- **License:** `custom <https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units>`__
- **Description:** Europe's NUTS administrative regions.

``data/bundle/nama_10r_3popgdp.tsv.gz``

- **Source:** Eurostat
- **Link:** https://ec.europa.eu/eurostat/databrowser/view/NAMA_10R_3POPGDP/default/table?lang=en
- **License:** `custom <https://ec.europa.eu/eurostat/about-us/policies/copyright>`__
- **Description:** Average annual population to calculate regional GDP data (thousand persons) by NUTS 3 regions.

``data/bundle/nama_10r_3gdp.tsv.gz``

- **Source:** Eurostat
- **Link:** https://ec.europa.eu/eurostat/databrowser/view/nama_10r_3gdp/default/table?lang=en
- **License:** `custom <https://ec.europa.eu/eurostat/about-us/policies/copyright>`__
- **Description:** Gross domestic product (GDP) at current market prices by NUTS 3 regions.

``data/bundle/corine``

- **Source:** European Environment Agency (EEA)
- **Link:** https://sdi.eea.europa.eu/catalogue/copernicus/api/records/a84ae124-c5c5-4577-8e10-511bfe55cc0d
- **License:** `custom <https://sdi.eea.europa.eu/catalogue/copernicus/api/records/a84ae124-c5c5-4577-8e10-511bfe55cc0d>`__
- **Description:** CORINE Land Cover (CLC) database.

``data/bundle/eea``

- **Source:** European Environment Agency (EEA)
- **Link:** https://www.eea.europa.eu/en/datahub/datahubitem-view/3b7fe76c-524a-439a-bfd2-a6e4046302a2
- **License:** CC-BY 4.0 (`reference <https://www.eea.europa.eu/en/legal-notice#copyright-notice>`__)
- **Description:** Total GHG emissions and removals in the EU.

``data/bundle/emobility``

- **Source:** German Federal Highway Research Institute (BASt)
- **Link:** https://www.bast.de/DE/Verkehrstechnik/Fachthemen/v2-verkehrszaehlung/zaehl_node.html
- **License:** CC-BY 4.0 (`reference <https://www.bast.de/DE/Verkehrstechnik/Fachthemen/v2-verkehrszaehlung/Nutzungsbedingungen.html?nn=1819490>`__)
- **Description:** Contains data from permanent automatic counting stations on highways and federal roads in Germany.

``data/bundle/h2_salt_caverns_GWh_per_sqkm.geojson``

- **Source:** Dilara Gulcin Caglayan, Nikolaus Weber, Heidi U. Heinrichs, Jochen
  Linßen, Martin Robinius, Peter A. Kukla, Detlef Stolten, Technical potential
  of salt caverns for hydrogen storage in Europe, International Journal of
  Hydrogen Energy, Volume 45, Issue 11, 2020, Pages 6793-6805.
- **Link:** https://doi.org/10.1016/j.ijhydene.2019.12.161
- **License:** CC-BY 4.0
- **Description:** Contains geological hydrogen storage potentials for Europe.

``data/bundle/natura``

- **Source:** European Environment Agency (EEA)
- **Link:** https://www.eea.europa.eu/en/datahub/datahubitem-view/6fc8ad2d-195d-40f4-bdec-576e7d1268e4
- **License:** CC-BY 4.0 (`reference <https://www.eea.europa.eu/en/legal-notice#copyright-notice>`__)
- **Description:** Natura 2000 natural protection areas.

``data/bundle/gebco``

- **Source:** GEBCO
- **Link:** https://www.gebco.net/data_and_products/gridded_bathymetry_data/version_20141103/
- **License:** CC0 (`reference <https://www.bodc.ac.uk/data/documents/nodb/301801/>`__)
- **Description:** Bathymetric dataset (2014 version).

``data/bundle/GDP_per_capita_PPP_1990_2015_v2.nc``

- **Source:** Kummu, M., Taka, M. & Guillaume, J. Gridded global datasets for
  Gross Domestic Product and Human Development Index over 1990-2015. Sci Data 5,
  180004 (2018). https://doi.org/10.1038/sdata.2018.4
- **Link:** https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0
- **License:** CC0 (`reference <https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0>`__)
- **Description:** Gridded GDP data.

``data/bundle/ppp_2013_1km_Aggregated.tif``

- **Source:** WorldPop (www.worldpop.org - School of Geography and Environmental
  Science, University of Southampton; Department of Geography and Geosciences,
  University of Louisville; Departement de Geographie, Universite de Namur) and
  Center for International Earth Science Information Network (CIESIN), Columbia
  University (2018). Global High Resolution Population Denominators Project -
  Funded by The Bill and Melinda Gates Foundation (OPP1134076).
  https://dx.doi.org/10.5258/SOTON/WP00647
- **Link:** https://hub.worldpop.org/doi/10.5258/SOTON/WP00647
- **License:** CC-BY 4.0 (`reference <https://hub.worldpop.org/geodata/summary?id=24770>`__)
- **Description:** Gridded population data.


Specific retrieval rules
========================

Data in this section is retrieved and extracted in rules specified in ``rules/retrieve.smk``.

``data/nuts``

- **Source:** GISCO
- **Link:** https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/
- **License:** `custom <https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units>`__
- **Description:** Europe's NUTS administrative regions.

``data/ENSPRESO_BIOMASS.xlsx``

- **Source:** European Commission Joint Research Centre (JRC)
- **Link:** https://data.jrc.ec.europa.eu/dataset/74ed5a04-7d74-4807-9eab-b94774309d9f
- **License:** CC-BY 4.0
- **Description:** Contains biomass potentials for Europe.

``data/complete_map_2020_unit_Mt.geojson``

- **Source:** SETIS
- **Link:** https://setis.ec.europa.eu/european-co2-storage-database_en, processed with https://github.com/ericzhou571/Co2Storage
- **License:** `various <https://setis.ec.europa.eu/european-co2-storage-database_en>`__
- **Description:** European CO2 storage database CO2StoP.

``data/myb1-2022-nitro-ert.xlsx``

- **Source:** United States Geological Survey (USGS)
- **Link:** https://www.usgs.gov/centers/national-minerals-information-center/nitrogen-statistics-and-information
- **License:** CC0 (`reference <https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits>`__)
- **Description:** Statistics and information on the worldwide supply of, demand for, and flow of the mineral commodity nitrogen.

``data/Industrial_Database.csv``

- **Source:** Simon Pezzutto, Stefano Zambotti, Silvia Croce, Pietro Zambelli,
  Giulia Garegnani, Chiara Scaramuzzino, Ramón Pascual Pascuas, Alyona
  Zubaryeva, Franziska Haas, Dagmar Exner (EURAC), Andreas Mueller (e-think),
  Michael Hartner (TUW), Tobias Fleiter, Anna-Lena Klingler, Matthias Kuehnbach,
  Pia Manz, Simon Marwitz, Matthias Rehfeldt, Jan Steinbach, Eftim Popovski
  (Fraunhofer ISI) Reviewed by Lukas Kranzl, Sara Fritz (TUW)
  Hotmaps Project, D2.3 WP2 Report - Open Data Set for the EU28, 2018
  https://www.hotmaps-project.eu
- **Link:** https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database
- **License:** CC-BY 4.0 (`reference <https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database>`__)
- **Description:** Contains georeferenced industrial sites of energy-intensive
  industry sectors, together with GHG-emissions, production capacity, fuel
  demand and excess heat potentials calculated from emission and production
  data.

``data/eurostat/Balances-April2023``

- **Source:** Eurostat
- **Link:** https://ec.europa.eu/eurostat/documents/38154/4956218/Balances-April2023.zip
- **License:** CC-BY 4.0 (`reference <https://commission.europa.eu/legal-notice_en>`__)
- **Description:** Contains energy balances for Europe.

``data/eurostat/eurostat-household_energy_balances-february_2024.csv``

- **Source:** Eurostat
- **Link:** https://ec.europa.eu/eurostat/databrowser-backend/api/extraction/1.0/LIVE/false/sdmx/csv/nrg_d_hhq__custom_11480365?startPeriod=2013&endPeriod=2022&i&compressed=true
- **License:** CC-BY 4.0 (`reference <https://commission.europa.eu/legal-notice_en>`__)
- **Description:** Contains household energy balances for Europe.

``data/jrc-idees-2021``

- **Source:** Rózsai, M., Jaxa-Rozen, M., Salvucci, R., Sikora, P., Tattini, J.
  and Neuwahl, F., JRC-IDEES-2021: the Integrated Database of the European
  Energy System - Data update and technical documentation, Publications Office
  of the European Union, Luxembourg, 2024, doi:10.2760/614599, JRC137809.
- **Link:** https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/JRC-IDEES/JRC-IDEES-2021_v1
- **License:** CC-BY 4.0 (`reference <https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/JRC-IDEES/copyright.txt>`__)
- **Description:** Contains more granular energy balances for Europe.

``data/gas_network``

- **Source:** Jan Diettrich, Adam Pluta, & Wided Medjroubi. (2021). SciGRID_gas
  IGGIELGN (1.1.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4767098
- **Link:** https://zenodo.org/records/4767098
- **License:** CC-BY 4.0 (`reference <https://zenodo.org/record/4767098>`__)
- **Description:** Contains gas infrastructure data.

``data/electricity_demand_raw.csv``

- **Source:** Open Power System Data (OPSD) from ENTSO-E Transparency
- **Link:**
  https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv
  and https://data.open-power-system-data.org/time_series/2020-10-06/time_series_60min_singleindex.csv
- **License:** unknown
- **Description:** Contains country-level electricity demand time series.

``data/load_synthetic_raw.csv``

- **Source:** Frysztacki, M., van der Most, L., & Neumann, F. (2024).
  Interannual Electricity Demand Calculator [Data set]. Zenodo.
  https://doi.org/10.5281/zenodo.10820928
- **Link:** https://zenodo.org/records/10820928
- **License:** CC-BY 4.0
- **Description:** Contains synthetic country-level electricity demand time series.

``data/shipdensity_global.zip``

- **Source:** World Bank
- **Link:** https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density
- **License:** CC-BY 4.0 (`reference <https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`__)
- **Description:** Global shipping traffic density.

``data/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif``

- **Source:** Marcel Buchhorn, Bruno Smets, Luc Bertels, Bert De Roo, Myroslava
  Lesiv, Nandin-Erdene Tsendbazar, Martin Herold, & Steffen Fritz. (2020).
  Copernicus Global Land Service: Land Cover 100m: collection 3: epoch 2019:
  Globe (V3.0.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3939050
- **Link:** https://zenodo.org/records/3939050
- **License:** CC-BY 4.0 (`reference <https://zenodo.org/record/3939050>`__)
- **Description:** Contains rastered land cover and land use data.

``data/LUISA_basemap_020321_50m.tif``

- **Source:** European Commission Joint Research Centre (JRC)
- **Link:** https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/LUISA/EUROPE/Basemaps/LandUse/2018/LATEST/
- **License:** CC-BY 4.0 (`reference <https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/LUISA/EUROPE/Basemaps/LandUse/2018/LATEST/>`__)
- **Description:** Contains rastered land cover and land use data.

``data/eez``

- **Source:** Marine Regions
- **Link:** https://www.marineregions.org/download_file.php
- **License:** CC-BY-NC-SA
- **Description:** Contains offshore exclusive economic zones.

``data/worldbank``

- **Source:** World Bank
- **Link:** https://data.worldbank.org/indicator/SP.URB.TOTL.IN.ZS
- **License:** CC-BY 4.0
- **Description:** Contains share of urban population by country.

``data/naturalearth``

- **Source:** Natural Earth
- **Link:** https://www.naturalearthdata.com/downloads/10m-cultural-vectors/
- **License:** CC0 (`reference <https://www.naturalearthdata.com/about/terms-of-use/>`__)
- **Description:** Country shapes, using point-of-view (POV) variant of Germany so that Crimea is included.

``data/gem/Europe-Gas-Tracker-2024-05.xlsx``

- **Source:** Global Energy Monitor
- **Link:** https://globalenergymonitor.org/projects/global-steel-plant-tracker/
- **License:** CC-BY 4.0 (`reference <https://globalenergymonitor.org/projects/europe-gas-tracker/download-data/>`__)
- **Description:** Covers methane gas pipelines, LNG terminals, oil and gas-fired power plants, and methane gas extraction sites.

``data/gem/Global-Steel-Plant-Tracker-April-2024-Standard-Copy-V1.xlsx``

- **Source:** Global Energy Monitor
- **Link:** https://globalenergymonitor.org/projects/global-steel-plant-tracker/
- **License:** CC-BY 4.0 (`reference <https://globalenergymonitor.org/projects/global-steel-plant-tracker/download-data/>`__)
- **Description:** The Global Steel Plant Tracker (GSPT) provides information on
  global crude iron and steel production plants, and includes every plant
  currently operating with a capacity of five hundred thousand tonnes per year
  (ttpa) or more of crude iron or steel.

``data/WDPA.gpkg``

- **Source:** UNEP-WCMC and IUCN (2024), Protected Planet: The World Database on
  Protected Areas (WDPA) [Online], September 2024, Cambridge, UK: UNEP-WCMC and
  IUCN. Available at: www.protectedplanet.net.
- **Link:** https://www.protectedplanet.net/en/thematic-areas/wdpa
- **License:** `custom <https://www.protectedplanet.net/en/legal>`__
- **Description:** Contains global protected areas.

``data/WDPA_WDOECM_marine.gpkg``

- **Source:** UNEP-WCMC and IUCN (2024), Protected Planet: The World Database on
  Protected Areas (WDPA) and World Database on Other Effective Area-based
  Conservation Measures (WD-OECM) [Online], September 2024, Cambridge, UK:
  UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.
- **Link:** https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
- **License:** `custom <https://www.protectedplanet.net/en/legal>`__
- **Description:** Contains global protected marine areas.

``data/osm-prebuilt``

- **Source:** OpenStreetMap; Xiong, B., Neumann, F., & Brown, T. (2024).
  Prebuilt Electricity Network for PyPSA-Eur based on OpenStreetMap Data (0.4)
  [Data set]. Zenodo. https://doi.org/10.5281/zenodo.13759222
- **Link:** https://zenodo.org/records/13759222
- **License:** ODbL (`reference <https://zenodo.org/records/13759222>`)
- **Description:** Pre-built data of high-voltage transmission grid in Europe from OpenStreetMap.

``data/osm-raw``

- **Source:** OpenStreetMap via Overpass API
- **Link:** https://overpass-api.de/api/interpreter
- **License:** ODbL
- **Description:** Data of high-voltage transmission grid in Europe from OpenStreetMap.

``cutouts``

- **Source:** `ERA5
  <https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview>`__
  and `SARAH-3 <https://navigator.eumetsat.int/product/EO:EUM:DAT:0863>`__
- **Link:** https://zenodo.org/records/12791128
- **License:** CC-BY 4.0
- **Description:** Contains weather data cutouts for Europe to read in with ``atlite``.

``resources/costs_{year}.csv``

- **Source:** various, mostly compiled from Danish Energy Agency (DEA)
  `Technology Catalogues
  <https://ens.dk/en/our-services/technology-catalogues>`__.
- **Link:** https://github.com/PyPSA/technology-data
- **License:** GPL-3.0
- **Description:** Contains technology data for different years such as costs, efficiencies, and lifetimes.

``resources/powerplants.csv``

- **Source:** F. Gotzens, H. Heinrichs, J. Hörsch, and F. Hofmann, Performing
  energy modelling exercises in a transparent way - The issue of data quality in
  power plant databases, Energy Strategy Reviews, vol. 23, pp. 1-12, Jan. 2019.
  https://doi.org/10.1016/j.esr.2018.11.004
- **Link:** https://github.com/PyPSA/powerplantmatching
- **License:** GPL-3.0
- **Description:** Contains matched dataset of powerplants in Europe.


Repository
==========

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
