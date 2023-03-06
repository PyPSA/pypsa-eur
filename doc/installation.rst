<<<<<<< HEAD
..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

=======
>>>>>>> pypsa-eur-sec/master
.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

<<<<<<< HEAD
Clone the Repository
====================

First of all, clone the `PyPSA-Eur repository <https://github.com/PyPSA/pypsa-eur>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces!
If you do not have ``git`` installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. code:: bash

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone https://github.com/PyPSA/pypsa-eur.git


.. _deps:

Install Python Dependencies
===============================

PyPSA-Eur relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to install them.
Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system.
For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are curated in the `envs/environment.yaml <https://github.com/PyPSA/pypsa-eur/blob/master/envs/environment.yaml>`_ file.
The environment can be installed and activated using
=======
Install PyPSA-Eur and its data
==============================

First install `PyPSA-Eur <https://github.com/PyPSA/pypsa-eur>`_ and all
its dependencies. Clone the repository:


.. code:: bash

    projects % git clone https://github.com/PyPSA/pypsa-eur.git

then download and unpack all the PyPSA-Eur data files by running the following snakemake rule:

.. code:: bash

    projects/pypsa-eur % snakemake -j 1 retrieve_databundle



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

You can create an environment using the environment.yaml file in pypsa-eur/envs:
>>>>>>> pypsa-eur-sec/master

.. code:: bash

    .../pypsa-eur % conda env create -f envs/environment.yaml

    .../pypsa-eur % conda activate pypsa-eur

<<<<<<< HEAD
Note that activation is local to the currently open shell!
After opening a new terminal window, one needs to reissue the second command!

.. note::
    If you have troubles with a slow ``conda`` installation, we recommend to install
    `mamba <https://github.com/QuantStack/mamba>`_ as a fast drop-in replacement via

    .. code:: bash

        conda install -c conda-forge mamba

    and then install the environment with

    .. code:: bash

        mamba env create -f envs/environment.yaml

Install a Solver
================

PyPSA passes the PyPSA-Eur network model to an external solver for performing a total annual system cost minimization with optimal power flow.
PyPSA is known to work with the free software

- `Ipopt <https://coin-or.github.io/Ipopt/INSTALL.html>`_
- `Cbc <https://projects.coin-or.org/Cbc#DownloadandInstall>`_
- `GLPK <https://www.gnu.org/software/glpk/>`_ (`WinGLKP <http://winglpk.sourceforge.net/>`_)
- `HiGHS <https://highs.dev/>`_

and the non-free, commercial software (for some of which free academic licenses are available)

- `Gurobi <https://www.gurobi.com/documentation/quickstart.html>`_
- `CPLEX <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_
- `FICO® Xpress Solver <https://www.fico.com/de/products/fico-xpress-solver>`_

For installation instructions of these solvers for your operating system, follow the links above.
Commercial solvers such as Gurobi and CPLEX currently significantly outperform open-source solvers for large-scale problems.
It might be the case that you can only retrieve solutions by using a commercial solver.

.. seealso::
    `Getting a solver in the PyPSA documentation <https://pypsa.readthedocs.io/en/latest/installation.html#getting-a-solver-for-linear-optimisation>`_

.. note::
    The rules :mod:`cluster_network` and :mod:`simplify_network` solve a quadratic optimisation problem for clustering.
    The open-source solvers Cbc and GlPK cannot handle this. A fallback to Ipopt is implemented in this case, but requires
    also Ipopt to be installed. For an open-source solver setup install in your ``conda`` environment on OSX/Linux

    .. code:: bash

        conda activate pypsa-eur
        conda install -c conda-forge ipopt coincbc

    and on Windows

    .. code:: bash

        conda activate pypsa-eur
        conda install -c conda-forge ipopt glpk

    or

    .. code:: bash

        conda activate pypsa-eur
        pip install highspy

.. _defaultconfig:

Set Up the Default Configuration
================================

PyPSA-Eur has several configuration options that must be specified in a ``config.yaml`` file located in the root directory.
An example configuration ``config.default.yaml`` is maintained in the repository.
More details on the configuration options are in :ref:`config`.

Before first use, create a ``config.yaml`` by copying the example.

.. code:: bash

    .../pypsa-eur % cp config.default.yaml config.yaml

Users are advised to regularly check their own ``config.yaml`` against changes in the ``config.default.yaml``
when pulling a new version from the remote repository.
=======
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
>>>>>>> pypsa-eur-sec/master
