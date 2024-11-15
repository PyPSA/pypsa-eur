..
  SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

#############
Base network
#############

.. raw:: html

   <iframe src="base-network-raw.html" width="100%" height="600px"></iframe>

The map might take a moment to load. To view it in full screen, click `here <base-network-raw.html>`__.

``data/osm-prebuilt``

- **Source:** OpenStreetMap; Xiong, B., Neumann, F., & Brown, T. (2024).
  Prebuilt Electricity Network for PyPSA-Eur based on OpenStreetMap Data (0.5)
  [Data set]. Zenodo. https://doi.org/10.5281/zenodo.13981528
- **Link:** https://zenodo.org/records/13981528
- **License:** ODbL (`reference <https://zenodo.org/records/13981528>`)
- **Description:** Pre-built data of high-voltage transmission grid in Europe from OpenStreetMap.

This dataset contains a topologically connected representation of the European
high-voltage grid (220 kV to 750 kV) constructed using OpenStreetMap data. Input data
was retrieved using the `Overpass turbo API <https://overpass-turbo.eu/>`__. A heurisitic
cleaning process was used to for lines and links where electrical parameters are
incomplete, missing, or ambiguous. Close substations within a radius of 500 m are
aggregated to single buses, exact locations of underlying substations is preserved.
Unique identifiers for lines and links are preserved.

A detailed explanation on the background, methodology, and validation of this dataset
can be found in `this paper <https://doi.org/10.48550/arXiv.2408.17178>`__ preprint
currently under peer-review.

Countries included in the dataset:

Albania (AL), Austria (AT), Belgium (BE), Bosnia and Herzegovina (BA), Bulgaria (BG),
Croatia (HR), Czech Republic (CZ), Denmark (DK), Estonia (EE), Finland (FI), France
(FR), Germany (DE), Greece (GR), Hungary (HU), Ireland (IE), Italy (IT), Kosovo (XK),
Latvia (LV), Lithuania (LT), Luxembourg (LU), Moldova (MD), Montenegro (ME), Netherlands
(NL), North Macedonia (MK), Norway (NO), Poland (PL), Portugal (PT), Romania (RO), Serbia
(RS), Slovakia (SK), Slovenia (SI), Spain (ES), Sweden (SE), Switzerland (CH), Ukraine
(UA), United Kingdom (GB)
