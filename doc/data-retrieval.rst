..
  SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

#########################
Specific retrieval rules
#########################

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
