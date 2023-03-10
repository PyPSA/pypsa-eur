..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##################
Cost Assumptions
##################

The database of cost assumptions is retrieved from the repository
`PyPSA/technology-data <https://github.com/pypsa/technology-data>`_ and then
saved to ``resources/costs.csv``.

The ``config.yaml`` provides options to choose a reference year (``costs: year:``) and use a specific version of the repository ``costs: version:``.

It includes cost assumptions for all included technologies for specific
years from various sources, namely for

- discount rate,
- lifetime,
- investment (CAPEX),
- fixed operation and maintenance (FOM),
- variable operation and maintenance (VOM),
- fuel costs,
- efficiency, and
- carbon-dioxide intensity.

The given overnight capital costs are annualised to net present costs
with a discount rate of :math:`r` over the economic lifetime :math:`n` using the annuity factor

.. math::

    a = \frac{1-(1+r)^{-n}}{r}.

Based on the parameters above the ``marginal_cost`` and ``capital_cost`` of the system components are calculated.


Modifying Cost Assumptions
==========================

Some cost assumptions (e.g. marginal cost and capital cost) can be directly overwritten in the ``config.yaml`` (cf. Section  :ref:`costs_cf`  in :ref:`config`).

To change cost assumptions in more detail, modify cost assumptions directly in ``resources/costs.csv`` as this is not yet supported through the config file.

You can also build multiple different cost databases. Make a renamed copy of ``resources/costs.csv`` (e.g. ``data/costs-optimistic.csv``) and set the variable ``COSTS=data/costs-optimistic.csv`` in the ``Snakefile``.
