.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

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

Preferred method: Pixi
----------------------

PyPSA-Eur relies on a set of other Python packages to function.
We manage these using [pixi](https://pixi.sh/latest/).
Once pixi is installed, you will have access to the PyPSA-Eur dependencies by prepending all your calls with `pixi run` (e.g. `pixi run snakemake -call -n`).
Alternatively, you can call `pixi shell` in your terminal to activate your working environment, then continue without needing to prepend calls with `pixi run` (e.g. `snakemake -call -n`).

Legacy method: conda
----------------------


If you cannot access `pixi` on your machine, you can also install using `conda`/`mamba`/`micromamba`.
To do so, you will install from one of our platform-specific environment files:

* For Intel/AMD processors:

  - Linux: ``envs/default_linux-64.pin.txt``
  - macOS: ``envs/default_osx-64.pin.txt``
  - Windows: ``envs/default_win-64.pin.txt``

* For ARM processors:

  - macOS (Apple Silicon): ``envs/default_osx-arm64.pin.txt``
  - Linux (ARM): Currently not supported via lock files; requires building certain packages, such as ``PySCIPOpt``, from source

.. code:: console

    $ conda update conda

    $ conda create -n pypsa-eur -f envs/default_linux-64.pin.txt # select the appropriate file for your platform

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
    The rules :mod:`cluster_network` solves a mixed-integer quadratic optimisation problem for clustering.
    The open-source solvers HiGHS, Cbc and GlPK cannot handle this. A fallback to SCIP is implemented in this case, which is included in the standard environment specifications.
    For an open-source solver setup install for example HiGHS **and** SCIP in your ``conda`` environment on OSX/Linux.
    To install the default solver Gurobi, run

    .. code:: console

        $ conda activate pypsa-eur
        $ conda install -c gurobi gurobi"=12.0.1"

    Additionally, you need to setup your `Gurobi license <https://www.gurobi.com/solutions/licensing/>`__.
