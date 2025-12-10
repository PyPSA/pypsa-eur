.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
Solving Networks
##########################################

After generating and clustering the networks, :mod:`compose_network` produces
``networks/composed_{horizon}.nc`` for each planning horizon. These files are
then passed to the single :mod:`solve_network` rule, which runs
``scripts/solve_network.py`` regardless of whether the study is electricity-only
or sector-coupled. Dispatch-only analyses on an already solved network remain
available through :mod:`solve_operations_network`.

.. _solve:

Rule ``solve_network``
=========================

.. automodule:: solve_network

.. _solve_operations:

Rule ``solve_operations_network``
====================================

.. automodule:: solve_operations_network

Rule ``make_summary``
=============================

.. automodule:: make_summary
