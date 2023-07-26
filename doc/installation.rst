..
  SPDX-FileCopyrightText: 2019-2023 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Clone the Repository
====================

First of all, clone the `PyPSA-Eur repository <https://github.com/PyPSA/pypsa-eur>`_ using the version control system ``git`` in the command line.

.. code:: bash

    /some/other/path % cd /some/path

    /some/path % git clone https://github.com/PyPSA/pypsa-eur.git


.. _deps:

Install Python Dependencies
===============================

PyPSA-Eur relies on a set of other Python packages to function.
We recommend using the package manager `mamba <https://mamba.readthedocs.io/en/latest/>`_ to install them and manage your environments.
For instructions for your operating system follow the ``mamba`` `installation guide <https://mamba.readthedocs.io/en/latest/installation.html>`_.
You can also use ``conda`` equivalently.

The package requirements are curated in the `envs/environment.yaml <https://github.com/PyPSA/pypsa-eur/blob/master/envs/environment.yaml>`_ file.
The environment can be installed and activated using

.. code:: bash

    .../pypsa-eur % mamba env create -f envs/environment.yaml

    .../pypsa-eur % mamba activate pypsa-eur

.. note::
    The equivalent commands for ``conda`` would be

    .. code:: bash

        .../pypsa-eur % conda env create -f envs/environment.yaml

        .../pypsa-eur % conda activate pypsa-eur


Install a Solver
================

PyPSA passes the PyPSA-Eur network model to an external solver for performing the optimisation.
PyPSA is known to work with the free software

- `HiGHS <https://highs.dev/>`_
- `Cbc <https://projects.coin-or.org/Cbc#DownloadandInstall>`_
- `GLPK <https://www.gnu.org/software/glpk/>`_ (`WinGLKP <http://winglpk.sourceforge.net/>`_)
- `Ipopt <https://coin-or.github.io/Ipopt/INSTALL.html>`_

and the non-free, commercial software (for some of which free academic licenses are available)

- `Gurobi <https://www.gurobi.com/documentation/quickstart.html>`_
- `CPLEX <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_
- `FICO Xpress Solver <https://www.fico.com/de/products/fico-xpress-solver>`_

For installation instructions of these solvers for your operating system, follow the links above.
Commercial solvers such as Gurobi and CPLEX currently significantly outperform open-source solvers for large-scale problems, and
it might be the case that you can only retrieve solutions by using a commercial solver.
Nevertheless, you can still use open-source solvers for smaller problems.

.. seealso::
    `Instructions how to install a solver in the documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/installation.html#getting-a-solver-for-linear-optimisation>`_

.. note::
    The rules :mod:`cluster_network` and :mod:`simplify_network` solve a quadratic optimisation problem for clustering.
    The open-source solvers Cbc and GlPK cannot handle this. A fallback to Ipopt is implemented in this case, but requires
    it to be installed. For an open-source solver setup install in your ``conda`` environment on OSX/Linux

    .. code:: bash

        mamba activate pypsa-eur
        mamba install -c conda-forge ipopt coincbc

    and on Windows

    .. code:: bash

        mamba activate pypsa-eur
        mamba install -c conda-forge ipopt glpk

    For HiGHS, run

    .. code:: bash

        mamba activate pypsa-eur
        mamba install -c conda-forge ipopt
        pip install highspy

    For Gurobi, run

    .. code:: bash

        mamba activate pypsa-eur
        mamba install -c gurobi gurobi

    Additionally, you need to setup your `Gurobi license <https://www.gurobi.com/solutions/licensing/>`_.


.. _defaultconfig:

Handling Configuration Files
============================

PyPSA-Eur has several configuration options that must be specified in a
``config/config.yaml`` file located in the root directory. An example configuration
``config/config.default.yaml`` is maintained in the repository, which will be used to
automatically create your customisable ``config/config.yaml`` on first use. More
details on the configuration options are in :ref:`config`.

You can also use ``snakemake`` to specify another file, e.g.
``config/config.mymodifications.yaml``, to update the settings of the ``config/config.yaml``.

.. code:: bash

    .../pypsa-eur % snakemake -call --configfile config/config.mymodifications.yaml

.. warning::
    Users are advised to regularly check their own ``config/config.yaml`` against changes
    in the ``config/config.default.yaml`` when pulling a new version from the remote
    repository.
