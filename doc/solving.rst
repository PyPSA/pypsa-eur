..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Solving Networks
##########################################

After generating and simplifying the networks they can be solved through the
rule :mod:`solve_network`  by using the collection rules ``solve_elec_networks``
or ``solve_sector_networks``. Moreover, networks can be solved for dispatch-only
analyses on an already solved network with :mod:`solve_operations_network`.

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

.. warning::
   More comprehensive documentation for this rule will be released soon.
