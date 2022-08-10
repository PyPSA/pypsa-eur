PyPSA-Eur-Sec: A Sector-Coupled Open Optimisation Model of the European Energy System
=====================================================================================

.. image:: https://img.shields.io/github/v/release/pypsa/pypsa-eur-sec?include_prereleases
    :alt: GitHub release (latest by date including pre-releases)

.. image:: https://readthedocs.org/projects/pypsa-eur-sec/badge/?version=latest
    :target: https://pypsa-eur-sec.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/github/license/pypsa/pypsa-eur-sec
    :alt: GitHub

.. image:: https://img.shields.io/github/repo-size/pypsa/pypsa-eur-sec
    :alt: GitHub repo size

.. image:: https://badges.gitter.im/PyPSA/community.svg
    :target: https://gitter.im/PyPSA/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
    :alt: Chat on Gitter


PyPSA-Eur-Sec is an open model dataset of the European energy system at the
transmission network level that covers the full ENTSO-E area.

PyPSA-Eur-Sec builds on the electricity generation and transmission
model `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ to add demand
and supply for the following sectors: transport, space and water
heating, biomass, energy consumption in the agriculture, industry 
and industrial feedstocks. This completes the energy system and includes 
all greenhouse gas emitters except waste management, agriculture, 
forestry and land use.


**WARNING**: PyPSA-Eur-Sec is under active development and has several
`limitations <https://pypsa-eur-sec.readthedocs.io/en/latest/limitations.html>`_ which
you should understand before using the model. The github repository
`issues <https://github.com/PyPSA/pypsa-eur-sec/issues>`_ collect known
topics we are working on (please feel free to help or make suggestions).
The `documentation <https://pypsa-eur-sec.readthedocs.io/>`_ remains somewhat
patchy.
We cannot support this model if you choose to use it.

.. note::
    You can find showcases of the model's capabilities in the
    preprint `Benefits of a Hydrogen Network in Europe
    <https://arxiv.org/abs/2207.05816>`_, a `paper in Joule with a
    description of the industry sector
    <https://arxiv.org/abs/2109.09563>`_, or in `a 2021 presentation
    at EMP-E <https://nworbmot.org/energy/brown-empe.pdf>`_.


This diagram gives an overview of the sectors and the links between
them:

.. image:: ../graphics/multisector_figure.png

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

Currently the scripts to solve and process the resulting PyPSA models
are also included in PyPSA-Eur-Sec, although they could in future be
better integrated with the corresponding scripts in PyPSA-Eur. A
stumbling block to sharing solve_network.py between PyPSA-Eur and
PyPSA-Eur-Sec is the different extra_functionality required to build
storage and CHP constraints.


PyPSA-Eur-Sec is designed to be imported into the open toolbox `PyPSA
<https://www.pypsa.org>`_ for which `documentation <https://pypsa.org/doc>`_ is
available as well.

This project is currently maintained by the `Department of Digital
Transformation in Energy Systems <https://tub-ensys.github.io>`_ at the
`Technical University of Berlin <https://www.tu.berlin>`_. Previous versions
were developed by the `Energy System Modelling group
<https://www.iai.kit.edu/english/2338.php>`_ at the `Institute for Automation
and Applied Informatics <https://www.iai.kit.edu/english/index.php>`_ at the
`Karlsruhe Institute of Technology <http://www.kit.edu/english/index.php>`_
which was funded by the `Helmholtz Association <https://www.helmholtz.de/en/>`_,
and by the `Renewable Energy Group
<https://fias.uni-frankfurt.de/physics/schramm/renewable-energy-system-and-network-analysis/>`_
at `FIAS <https://fias.uni-frankfurt.de/>`_ to carry out simulations for the
`CoNDyNet project <http://condynet.de/>`_, financed by the `German Federal
Ministry for Education and Research (BMBF) <https://www.bmbf.de/en/index.html>`_
as part of the `Stromnetze Research Initiative
<http://forschung-stromnetze.info/projekte/grundlagen-und-konzepte-fuer-effiziente-dezentrale-stromnetze/>`_.


Documentation
=============

**Getting Started**

* :doc:`installation`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Getting Started

   installation

**Implementation details**

* :doc:`spatial_resolution`
* :doc:`supply_demand`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Implementation details

   spatial_resolution
   supply_demand


**Foresight options**

* :doc:`overnight`
* :doc:`myopic`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Foresight options

   overnight
   myopic

**References**

* :doc:`release_notes`
* :doc:`limitations`

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: References

   release_notes
   limitations


Licence
=======

The code in PyPSA-Eur-Sec is released as free software under the
`MIT license <https://opensource.org/licenses/MIT>`_, see
`LICENSE <https://github.com/PyPSA/pypsa-eur-sec/blob/master/LICENSE.txt>`_.
However, different licenses and terms of use may apply to the various input data.
