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

Data in this section is downloaded and extracted from the Zenodo data bundle (https://zenodo.org/records/12760663).

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

``data/bundle/nuts``

- **Source:** GISCO
- **Link:** https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/
- **License:** `custom <https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units>`__
- **Description:** Europe's NUTS administrative regions.

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

``data/bunlde/ppp_2013_1km_Aggregated.tif``

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

``https://zenodo.org/records/10356004/files/ENSPRESO_BIOMASS.xlsx``

- **Source:** European Commission Joint Research Centre (JRC)
- **Link:** https://data.jrc.ec.europa.eu/dataset/74ed5a04-7d74-4807-9eab-b94774309d9f
- **License:** CC-BY 4.0
- **Description:** Contains biomass potentials for Europe.

``https://raw.githubusercontent.com/ericzhou571/Co2Storage/main/resources/complete_map_2020_unit_Mt.geojson``

- **Source:** SETIS
- **Link:** https://setis.ec.europa.eu/european-co2-storage-database_en, processed with https://github.com/ericzhou571/Co2Storage
- **License:** `various <https://setis.ec.europa.eu/european-co2-storage-database_en>`__
- **Description:** European CO2 storage database CO2StoP.

``https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/myb1-2022-nitro-ert.xlsx``

- **Source:** United States Geological Survey (USGS)
- **Link:** https://www.usgs.gov/centers/national-minerals-information-center/nitrogen-statistics-and-information
- **License:** CC0 (`reference <https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits>`__)
- **Description:** Statistics and information on the worldwide supply of, demand for, and flow of the mineral commodity nitrogen.

``https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database/-/raw/master/data/Industrial_Database.csv``

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

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/eurostat/eurostat-household_energy_balances-february_2024.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/jrc-idees-2021``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/gas_network``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/electricity_demand_raw.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/load_synthetic_raw.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/shipdensity_global.zip``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/LUISA_basemap_020321_50m.tif``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/eez``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/worldbank``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/naturalearth``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/gem/Europe-Gas-Tracker-2024-05.xlsx``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/gem/Global-Steel-Plant-Tracker-April-2024-Standard-Copy-V1.xlsx``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/WDPA.gpkg``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/WDPA_WDOECM_marine.gpkg``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/osm-prebuilt``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/osm-raw``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``cutouts``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``resources/costs_{year}.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``resources/powerplants.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.


Repository
==========

Data in this section is included in the PyPSA-Eur repository in the ``data`` folder.

``data/entsoegridkit``

- **Source:** ENTSO-E
- **Link:** https://www.entsoe.eu/data/map/, extracted with https://github.com/PyPSA/GridKit/tree/master/entsoe
- **License:** unknown
- **Description:** Data of high-voltage transmission grid in Europe.

``data/existing_infrastructure``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/retro``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/transmission_projects``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

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

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/era5-annual-runoff-per-country.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/gr-e-11.03.02.01.01-cc.csv``

- **Source:** Swiss Federal Statistics Office
- **Link:** https://www.bfs.admin.ch/asset/de/30305426
- **License:** `custom (OPEN BY ASK) <https://www.bfs.admin.ch/bfs/en/home/fso/swiss-federal-statistical-office/terms-of-use.html>`__
- **Description:** Stock of road motor vehicles in Switzerland.

``data/heat_load_profile_BDEW.csv``

- **Source:** oemof/demandlib
- **Link:** https://github.com/oemof/demandlib
- **License:** MIT
- **Description:** Contains standard heat load profiles based on data from BDEW (German Association of Energy and Water Industries).

.. note::
   The provenance of the data is unclear. We will improve this in the future.

``data/hydro_capacities.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

.. note::
   The provenance of the data is unclear. We will improve this in the future.

``data/links_p_nom.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

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

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/unit_commitment.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.

``data/biomass_transport_costs_supply_chain{1,2}.csv``

- **Source:**
- **Link:** https://www.xyz.com/renewables
- **License:** CC-BY 4.0
- **Description:** Contains energy balances for Europe.
