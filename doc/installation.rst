.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Install PyPSA-Eur and its data
==============================

First install `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ and all
its dependencies. Clone the repository:


.. code:: bash

    projects % git clone https://github.com/PyPSA/pypsa-eur.git

then download and unpack all the PyPSA-Eur data files by running the following snakemake rule:

.. code:: bash

    projects/pypsa-eur % snakemake -j 1 retrieve_databundle


Clone technology-data repository
================================

Next install the technology assumptions database `technology-data <https://github.com/PyPSA/technology-data>`_ by creating a parallel directory:

.. code:: bash

    projects % git clone https://github.com/PyPSA/technology-data.git


Clone PyPSA-Eur-Sec repository
==============================

Create a parallel directory for `PyPSA-Eur-Sec <https://github.com/PyPSA/pypsa-eur-sec>`_ with:

.. code:: bash

    projects % git clone https://github.com/PyPSA/pypsa-eur-sec.git

Environment/package requirements
================================



The requirements are the same as `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_. For
``solve_network.py`` in addition you need ``gurobipy``.  If you have
xarray version >= 0.15.1, you will need the latest master branch of
atlite version 0.0.2.

You can create an enviroment using the environment.yaml file in pypsa-eur/envs:

.. code:: bash

    .../pypsa-eur % conda env create -f envs/environment.yaml

    .../pypsa-eur % conda activate pypsa-eur

See details in `PyPSA-Eur Installation <https://pypsa-eur.readthedocs.io/en/latest/installation.html>`_

Data requirements
=================

Small data files are included directly in the git repository, while
larger ones are archived in a data bundle on zenodo (`10.5281/zenodo.5824485 <https://doi.org/10.5281/zenodo.5824485>`_).
The data bundle's size is around 640 MB.

To download and extract the data bundle on the command line:

.. code:: bash

    projects/pypsa-eur-sec/data % wget "https://zenodo.org/record/5824485/files/pypsa-eur-sec-data-bundle.tar.gz"
    projects/pypsa-eur-sec/data % tar -xvzf pypsa-eur-sec-data-bundle.tar.gz


The data licences and sources are given in the following table.


.. csv-table::
   :header-rows: 1
   :file: data.csv



Set up the default configuration
================================

First make your own copy of the ``config.yaml`` based on 
 ``config.default.yaml``. For example:

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

    projects/pypsa-eur-sec % snakemake -j1

Warning: you may need a computer cluster for this (with e.g. 10-100 GB of RAM
and several processors).

To only prepare the networks, you can run the scripts up to the point before optimization:

.. code:: bash

    projects/pypsa-eur-sec % snakemake -j1 prepare_sector_networks
