..
  SPDX-FileCopyrightText: 2019-2024 The PyPSA-Spain Authors

  SPDX-License-Identifier: CC-BY-4.0


PyPSA-Spain: an extension of PyPSA-Eur to model the Spanish energy system
=======================================================================================

PyPSA-Spain is an open-source model of the Spanish energy system based on the European model PyPSA-Eur.
The primary motivation behind the development of PyPSA-Spain was to leverage the benefits of a national energy model over a regional one, like the availability of specific datasets from national organisations. Additionally, a single-country model enables higher spatial and temporal resolution with the same computational resources due to the smaller geographical domain. Finally, it does not require assumptions about coordinated action between countries, making it a more suitable tool for analysing national energy policies. To accommodate cross-border interactions, a nested model approach with PyPSA-Eur was used, wherein electricity prices from neighbouring countries are precomputed through the optimisation of the European energy system.

PyPSA-Spain is an up-to-date fork of PyPSA-Eur, ensuring that advancements and bug fixes made to PyPSA-Eur are integrated. In addition, PyPSA-Spain includes a number of novel functionalities that enhance the representation of the Spanish energy system, as compared with PyPSA-Eur. 

Find the details in `this preprint introducing PyPSA-Spain <https://arxiv.org/abs/2412.06571>`__.

Further details will be posted soon, stay tuned..



.. image:: img/base.jpg
    :width: 70%
    :align: center



.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Getting Started

   introduction
   configuration
