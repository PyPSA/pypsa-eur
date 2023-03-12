This directory contains the TYNDP2020 which was extended and modified such that it can be easily integrated into PyPSA-Eur:

- Added bus coordinates
- Split assets into buses, lines and links
- Classified assets as 'new' or 'upgraded' with respect to the latest PyPSA-Eur gridextract from 08-03-2022 by performing coordinate-based matching of TYNDP-assets to the gridextract assets (approach taken from integration of TYNDP2018 links in `scripts/base_network.py`)
  - 'Upgraded' assets: Adopted bus coordinates and IDs from gridextract, took larger values of significant fields (e.g. voltage) to ensure that nothing is downgraded.
- Convert data to gridextract format

For more details on the extension, modification and merging see the external project [tyndp_to_pypsa](https://github.com/grecht/tyndp_to_pypsa) where this was implemented.
If the integration into PyPSA-Eur fails this could be due to changes in the gridextract or the file ´data/links_tyndp.csv´.
In that case update this data in  [tyndp_to_pypsa](https://github.com/grecht/tyndp_to_pypsa) and follow the instructions there to create the files in this directory again.
