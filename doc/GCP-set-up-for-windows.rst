..
  SPDX-FileCopyrightText: 2020 Maximilian Parzen and Emmanuel Paez
  
  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Set-up Google Cloud Platform on Windows (300$ free trial)
##########################################

Purpose
====================
The Google Cloud Platform is an ideal tool to test PyPSA-Eur or PyPSA-Eur-Sec when, 

- you don't have immediately access to a high performance computation facility,
- you have problems with the Windows operating system and want a quick run on a linux-based system,
- you want to model whole-EU in acceptable spatial and time-resolution to run small research or work projects,
- you need quick results (the GCP provide you in the trial version max. 32 vCPU cores and more than 600 GB memory)

What you basically do with the Google Cloud Platform is that you set up a virtual machine/computer in the cloud which can store and operate data.
Similar as on your local computer, you have to install all software and solvers, and create paths on the virtual machine to run PyPSA. 
The 300$ free trial google budget for the first Google Cloud Platform use equals roughly, 10 whole-EU simulations, with 181 nodes at hourly basis.

The following steps are required for a successfull Google CLoud Platform set-up:

- `Google Cloud Platform registration <https://console.cloud.google.com>`_, to receive 300$ free budget.
- `Creating an Virtual Machine (VM) instance <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_, which is practically a virtual computer with Linux as OS.
- `Installation of Cloud SDK <https://cloud.google.com/sdk/>`_, to create a communication channel between your computer and the cloud virtual machine (VM).
- `Installation of WinSCP <https://winscp.net/eng/download.php>`_, to comfortably handle or transfer files between the VM and you local computer.

- `Installation of PuTTy <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_ ## not sure about that one


Clone the Repository
====================

First of all, clone the `PyPSA-Eur repository <https://github.com/PyPSA/pypsa-eur>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces!

.. code:: bash

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone https://github.com/PyPSA/pypsa-eur.git

.. note::
    If you do not have ``git`` installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. _deps:

Install Python Dependencies
===============================

PyPSA-Eur relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to install them.
Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system.
For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are curated in the `environment.yaml <https://github.com/PyPSA/pypsa-eur/blob/master/environment.yaml>`_ file.
The environment can be installed and activated using

.. code:: bash

    .../pypsa-eur % conda env create -f environment.yaml

    .../pypsa-eur % conda activate pypsa-eur

.. note::
    Note that activation is local to the currently open shell!
    After opening a new terminal window, one needs to reissue the second command!

.. note::
    If you have troubles with a slow ``conda`` installation, we recommend to install
    `mamba <https://github.com/QuantStack/mamba>`_ as a fast drop-in replacement via

    .. code:: bash
        
        conda install -c conda-forge mamba

    and then install the environment with

    .. code:: bash

        mamba env create -f environment.yaml

Install a Solver
================

PyPSA passes the PyPSA-Eur network model to an external solver for performing a total annual system cost minimization with optimal power flow.
PyPSA is known to work with the free software

- `Ipopt <https://coin-or.github.io/Ipopt/INSTALL.html>`_
- `Cbc <https://projects.coin-or.org/Cbc#DownloadandInstall>`_
- `GLPK <https://www.gnu.org/software/glpk/>`_ (`WinGLKP <http://winglpk.sourceforge.net/>`_)

and the non-free, commercial software (for which free academic licenses are available)

- `Gurobi <https://www.gurobi.com/documentation/quickstart.html>`_
- `CPLEX <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_

and any other solver that works with the underlying modelling framework `Pyomo <http://www.pyomo.org/>`_.
For installation instructions of these solvers for your operating system, follow the links above.

.. seealso::
    `Getting a solver in the PyPSA documentation <https://pypsa.readthedocs.io/en/latest/installation.html#getting-a-solver-for-linear-optimisation>`_

.. note::
    Commercial solvers such as Gurobi and CPLEX currently significantly outperform open-source solvers for large-scale problems.
    It might be the case that you can only retrieve solutions by using a commercial solver.

.. note::
    The rules :mod:`cluster_network` and :mod:`simplify_network` solve a quadratic optimisation problem for clustering.
    The open-source solvers Cbc and GlPK cannot handle this. A fallback to Ipopt is implemented in this case, but requires
    also Ipopt to be installed. For an open-source solver setup install in your `conda` environment on OSX/Linux

    .. code:: bash

        conda activate pypsa-eur
        conda install -c conda-forge ipopt coincbc

    and on Windows

    .. code:: bash
        
        conda activate pypsa-eur
        conda install -c conda-forge ipopt glpk
        

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

.. Using PyPSA-Eur with Docker Images
.. ==================================

.. If docker. Optional.
.. To run on cloud computing.
.. Gurobi license - floating token server - license must not be tied to a particular machine
.. Provide ``Dockerfile``.
