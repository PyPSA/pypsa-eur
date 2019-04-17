# PyPSA-Eur-Sec: A Sector-Coupled Open Optimisation Model of the European Energy System

PyPSA-Eur-Sec builds on the electricity generation and transmission
model [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) to add demand
and supply for the following sectors: transport, space and water
heating, biomass, industry and industrial feedstocks. This completes
the energy system and includes all greenhouse gas emitters except
waste management, agriculture, forestry and land use.

PyPSA-Eur-Sec includes PyPSA-Eur as a
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
[subworkflow](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-sub-workflows). PyPSA-Eur-Sec
uses PyPSA-Eur to build the clustered transmission model along with
wind, solar PV and hydroelectricity potentials and time series. Then
PyPSA-Eur-Sec adds other conventional generators, storage units and
the additional sectors.

Currently the scripts to solve and process the resulting PyPSA models
are also included in PyPSA-Eur-Sec, although they could in future be
better integrated with the corresponding scripts in PyPSA-Eur.

# Installation

First install [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) and all
its dependencies.

Create a parallel directory for PyPSA-Eur-Sec with:
```shell
projects % git clone git@github.com:nworbmot/pypsa-eur-sec.git
```



## Data requirements

The Data requirements include JRC-IDEES-2015, JRC biomass potentials,
EEA emission statistics, Eurostat Energy Balances, urban district
heating potentials, emobility statistics, timezone mappings and
heating profiles.

```shell
projects/pypsa-eur-sec/data % wget "https://nworbmot.org/pypsa-eur-sec-data-bundle-190417.tar.gz"
projects/pypsa-eur-sec/data % tar xvzf pypsa-eur-sec-data-bundle-190417.tar.gz
```
