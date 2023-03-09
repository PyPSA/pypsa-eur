..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Solving Networks
##########################################

After generating and simplifying the networks they can be solved through the rule :mod:`solve_network`  by using the collection rule :mod:`solve_all_networks`. Moreover, networks can be solved for another focus with the derivative rules :mod:`solve_network`  by using the collection rule :mod:`solve_operations_network` for dispatch-only analyses on an already solved network.

.. _solve:

Rule ``solve_network``
=========================

.. automodule:: solve_network

.. _solve_operations:

Rule ``solve_operations_network``
====================================

.. automodule:: solve_operations_network

Rule ``solve_sector_network``
=============================

.. automodule:: solve_sector_network