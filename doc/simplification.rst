

..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Simplifying Electricity Networks
##########################################

The simplification ``snakemake`` rules prepare **approximations** of the full model, for which it is computationally viable to co-optimize generation, storage and transmission capacities.

- :mod:`simplify_network` transforms the transmission grid to a 380 kV only equivalent network, while
- :mod:`cluster_network` uses a `k-means <https://en.wikipedia.org/wiki/K-means_clustering>`_ based clustering technique to partition the network into a given number of zones and then reduce the network to a representation with one bus per zone.

The simplification and clustering steps are described in detail in the paper

- Jonas Hörsch and Tom Brown. `The role of spatial scale in joint optimisations of generation and transmission for European highly renewable scenarios <https://arxiv.org/abs/1705.07617>`_), *14th International Conference on the European Energy Market*, 2017. `arXiv:1705.07617 <https://arxiv.org/abs/1705.07617>`_, `doi:10.1109/EEM.2017.7982024 <https://doi.org/10.1109/EEM.2017.7982024>`_.

After simplification and clustering of the network, additional components may be appended in the rule :mod:`add_extra_components` and the network is prepared for solving in :mod:`prepare_network`.

.. _simplify:

Rule ``simplify_network``
============================

.. automodule:: simplify_network

.. _cluster:

Rule ``cluster_network``
===========================

.. automodule:: cluster_network

.. _extra_components:

Rule ``add_extra_components``
=============================

.. automodule:: add_extra_components


.. _prepare:

Rule ``prepare_network``
===========================

.. automodule:: prepare_network
