![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/pypsa/pypsa-eur-sec?include_prereleases)
[![Documentation](https://readthedocs.org/projects/pypsa-eur-sec/badge/?version=latest)](https://pypsa-eur-sec.readthedocs.io/en/latest/?badge=latest)
![GitHub](https://img.shields.io/github/license/pypsa/pypsa-eur-sec)
![Size](https://img.shields.io/github/repo-size/pypsa/pypsa-eur-sec)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3938042.svg)](https://doi.org/10.5281/zenodo.3938042)
[![Gitter](https://badges.gitter.im/PyPSA/community.svg)](https://gitter.im/PyPSA/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

# PyPSA-Eur-Sec: A Sector-Coupled Open Optimisation Model of the European Energy System

PyPSA-Eur-Sec is an open model dataset of the European energy system at the
transmission network level that covers the full ENTSO-E area.

PyPSA-Eur-Sec builds on the electricity generation and transmission
model [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) to add demand
and supply for the following sectors: transport, space and water
heating, biomass, industry and industrial feedstocks, agriculture,
forestry and fishing. This completes the energy system and includes
all greenhouse gas emitters except waste management and land use.

**WARNING**: PyPSA-Eur-Sec is under active development and has several
[limitations](https://pypsa-eur-sec.readthedocs.io/en/latest/limitations.html) which
you should understand before using the model. The github repository
[issues](https://github.com/PyPSA/pypsa-eur-sec/issues) collect known
topics we are working on (please feel free to help or make suggestions).
The [documentation](https://pypsa-eur-sec.readthedocs.io/) remains somewhat
patchy.
You can find showcases of the model's capabilities in the preprint
[Benefits of a Hydrogen Network in Europe](https://arxiv.org/abs/2207.05816),
a [paper in Joule with a description of the industry
sector](https://arxiv.org/abs/2109.09563), or in [a 2021
presentation at EMP-E](https://nworbmot.org/energy/brown-empe.pdf).
We cannot support this model if you choose to use it.

Please see the [documentation](https://pypsa-eur-sec.readthedocs.io/)
for installation instructions and other useful information about the snakemake workflow.

This diagram gives an overview of the sectors and the links between
them:

![sector diagram](graphics/multisector_figure.png)

Each of these sectors is built up on the transmission network nodes
from [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur):

![network diagram](https://github.com/PyPSA/pypsa-eur/blob/master/doc/img/base.png?raw=true)

For computational reasons the model is usually clustered down
to 50-200 nodes.


PyPSA-Eur-Sec was initially based on the model PyPSA-Eur-Sec-30 described
in the paper [Synergies of sector coupling and transmission
reinforcement in a cost-optimised, highly renewable European energy
system](https://arxiv.org/abs/1801.05290) (2018) but it differs by
being based on the higher resolution electricity transmission model
[PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) rather than a
one-node-per-country model, and by including biomass, industry,
industrial feedstocks, aviation, shipping, better carbon management,
carbon capture and usage/sequestration, and gas networks.


PyPSA-Eur-Sec includes PyPSA-Eur as a
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
[subworkflow](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-sub-workflows). PyPSA-Eur-Sec
uses PyPSA-Eur to build the clustered transmission model along with
wind, solar PV and hydroelectricity potentials and time series. Then
PyPSA-Eur-Sec adds other conventional generators, storage units and
the additional sectors.


# Licence

The code in PyPSA-Eur-Sec is released as free software under the
[MIT License](https://opensource.org/licenses/MIT), see `LICENSE.txt`.
However, different licenses and terms of use may apply to the various
input data.
