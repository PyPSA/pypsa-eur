.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Install PyPSA-Eur
=================

First install `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ and all
its dependencies. Clone the repository:


.. code:: bash

    projects % git clone git@github.com:PyPSA/pypsa-eur.git

then download and unpack all the PyPSA-Eur data files.


Clone PyPSA-Eur-Sec repository
==============================

Create a parallel directory for PyPSA-Eur-Sec with:

.. code:: bash

    projects % git clone git@github.com:PyPSA/pypsa-eur-sec.git

Environment/package requirements
================================



The requirements are the same as `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_. For
``solve_network.py`` in addition you need ``gurobipy``.  If you have
xarray version >= 0.15.1, you will need the latest master branch of
atlite version 0.0.2.


Data requirements
=================

The data requirements include the JRC-IDEES-2015 database, JRC biomass
potentials, EEA emission statistics, Eurostat Energy Balances, urban
district heating potentials, emobility statistics, timezone mappings
and heating profiles.

The data bundle is about 640 MB.

To download and extract it on the command line:

.. code:: bash

    projects/pypsa-eur-sec/data % wget "https://nworbmot.org/pypsa-eur-sec-data-bundle-190719.tar.gz"
    projects/pypsa-eur-sec/data % tar xvzf pypsa-eur-sec-data-bundle-190719.tar.gz

Set up the default configuration
================================

First make your own copy of the ``config.yaml``:

.. code:: bash

    projects/pypsa-eur-sec % cp config.default.yaml config.yaml


Getting started
===============


In ``config.yaml`` you can control the settings for the scenarios you
want to run, such as the number of nodes, the CO2 limit, the
installable potentials for solar and wind, which technologies are
activated, etc.

To run the full optimization with your settings:

.. code:: bash

    projects/pypsa-eur-sec % snakemake

Warning: you may need a computer cluster for this (with e.g. 10-100 GB of RAM
and several processors).

To only prepare the networks, you can run the scripts up to the point before optimization:

.. code:: bash

    projects/pypsa-eur-sec % snakemake prepare_sector_networks
