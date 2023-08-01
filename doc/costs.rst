..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

############################
Techno-Economic Assumptions
############################

The database of cost assumptions is retrieved from the repository
`PyPSA/technology-data <https://github.com/pypsa/technology-data>`_ and then
saved to a file ``data/costs_{year}.csv``. The ``config/config.yaml`` provides options
to choose a reference year and use a specific version of the repository.

.. literalinclude:: ../config/config.default.yaml
   :language: yaml
   :start-at: costs:
   :end-at:   version:

The file includes cost assumptions for all included technologies for specific
years compiled from various sources, namely for

- discount rate,
- lifetime,
- investment (CAPEX),
- fixed operation and maintenance (FOM),
- variable operation and maintenance (VOM),
- fuel costs,
- efficiency, and
- carbon-dioxide intensity.

Many values are taken from a database published by the Danish Energy Agency (`DEA
<https://ens.dk/en/our-services/projections-and-models/technology-data>`_).


The given overnight capital costs are annualised to net present costs
with a discount rate of :math:`r` over the economic lifetime :math:`n` using the annuity factor

.. math::

    a = \frac{1-(1+r)^{-n}}{r}.

Based on the parameters above the ``marginal_cost`` and ``capital_cost`` of the
system components are automatically calculated.


Modifying Assumptions
=====================

Some cost assumptions (e.g. marginal cost and capital cost) can be directly
set in the ``config/config.yaml`` (cf. Section  :ref:`costs_cf`  in
:ref:`config`). To change cost assumptions in more detail, make a copy of
``data/costs_{year}.csv`` and reference the new cost file in the ``Snakefile``:

.. literalinclude:: ../Snakefile
   :start-at: COSTS
   :end-at: COSTS
