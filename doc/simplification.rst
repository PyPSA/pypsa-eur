

..
  SPDX-FileCopyrightText: 2019-2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Simplifying Electricity Networks
##########################################

The simplification ``snakemake`` rules prepare **approximations** of the network model, for which it is computationally viable to co-optimize generation, storage and transmission capacities.

- :mod:`simplify_network` transforms the transmission grid to a 380 kV only equivalent network, while
- :mod:`cluster_network` uses a `k-means <https://en.wikipedia.org/wiki/K-means_clustering>`__ based clustering technique to partition the network into a given number of zones and then reduce the network to a representation with one bus per zone.

The simplification and clustering steps are described in detail in the paper

- Jonas Hörsch and Tom Brown. `The role of spatial scale in joint optimisations of generation and transmission for European highly renewable scenarios <https://arxiv.org/abs/1705.07617>`__), *14th International Conference on the European Energy Market*, 2017. `arXiv:1705.07617 <https://arxiv.org/abs/1705.07617>`__, `doi:10.1109/EEM.2017.7982024 <https://doi.org/10.1109/EEM.2017.7982024>`__.

After simplification and clustering of the network, further electricity network components may be appended in the rule :mod:`add_electricity` and the network is prepared for solving in :mod:`prepare_network`.

.. _simplify:

Rule ``simplify_network``
============================

.. automodule:: simplify_network

.. _cluster:

Rule ``cluster_network``
===========================

.. automodule:: cluster_network

.. _electricity:

Rule ``add_electricity``
=============================

.. automodule:: add_electricity

.. _prepare:

Rule ``prepare_network``
===========================

.. automodule:: prepare_network
