..
  SPDX-FileCopyrightText: 2019-2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Clone the Repository
====================

First of all, clone the `PyPSA-Eur repository <https://github.com/PyPSA/pypsa-eur>`__ using the version control system ``git`` in the command line.

.. code:: console

    $ git clone https://github.com/PyPSA/pypsa-eur.git


.. _deps:

Install Python Dependencies
===============================

PyPSA-Eur relies on a set of other Python packages to function.
We recommend using the package manager `mamba <https://mamba.readthedocs.io/en/latest/>`__
to install them and manage your environments. For instructions for your operating
system follow the ``mamba`` `installation guide <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`__.
You can also use ``conda`` equivalently.

The package requirements are curated in the ``envs/environment.yaml`` file.
There are also regularly updated pinned environment files for each operating system to
ensure reproducibility (``envs/windows-pinned.yaml``, ``envs/linux-pinned.yaml``, ``envs/macos-pinned.yaml``).
We recommend to use the pinned files for a stable environment, but you could also use
the unpinned file.

.. code:: console

    $ mamba env create -f envs/linux-pinned.yaml # replace for your os

    $ mamba activate pypsa-eur

.. note::
    The equivalent commands for ``conda`` would be

    .. code:: console

        $ conda env create -f envs/linux-pinned.yaml # replace for your os

        $ conda activate pypsa-eur


Install a Solver
================

PyPSA passes the PyPSA-Eur network model to an external solver for performing the optimisation.
PyPSA is known to work with the free software

- `HiGHS <https://highs.dev/>`__
- `Cbc <https://projects.coin-or.org/Cbc#DownloadandInstall>`__
- `GLPK <https://www.gnu.org/software/glpk/>`__ (`WinGLKP <http://winglpk.sourceforge.net/>`__)
- `SCIP <https://scipopt.github.io/PySCIPOpt/docs/html/index.html>`__

and the non-free, commercial software (for some of which free academic licenses are available)

- `Gurobi <https://www.gurobi.com/documentation/quickstart.html>`__
- `CPLEX <https://www.ibm.com/products/ilog-cplex-optimization-studio>`__
- `FICO Xpress Solver <https://www.fico.com/de/products/fico-xpress-solver>`__

For installation instructions of these solvers for your operating system, follow the links above.
Commercial solvers such as Gurobi and CPLEX currently significantly outperform open-source solvers for large-scale problems, and
it might be the case that you can only retrieve solutions by using a commercial solver.
Nevertheless, you can still use open-source solvers for smaller problems.

.. seealso::
    `Instructions how to install a solver in the documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/installation.html#getting-a-solver-for-linear-optimisation>`__

.. note::
    The rules :mod:`cluster_network` and :mod:`simplify_network` solve a mixed-integer quadratic optimisation problem for clustering.
    The open-source solvers HiGHS, Cbc and GlPK cannot handle this. A fallback to SCIP is implemented in this case, which is included in the standard environment specifications.
    For an open-source solver setup install for example HiGHS **and** SCIP in your ``conda`` environment on OSX/Linux.
    To install the default solver Gurobi, run

    .. code:: console

        $ mamba activate pypsa-eur
        $ mamba install -c gurobi gurobi

    Additionally, you need to setup your `Gurobi license <https://www.gurobi.com/solutions/licensing/>`__.


.. _defaultconfig:

Handling Configuration Files
============================

PyPSA-Eur has several configuration options that users can specify in a
``config/config.yaml`` file. The default configuration
``config/config.default.yaml`` is maintained in the repository. More details on
the configuration options are in :ref:`config`.

You can also use ``snakemake`` to specify another file, e.g.
``config/config.mymodifications.yaml``, to update the settings of the ``config/config.yaml``.

.. code:: console

    $ snakemake --configfile config/config.mymodifications.yaml
