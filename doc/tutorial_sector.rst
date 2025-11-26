.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _tutorial_sector:

###############################
Tutorial: Sector-Coupled
###############################

.. note::
    If you have not done it yet, follow the :ref:`installation` steps first.

    Also, checkout the tutorial for electricity-only systems first at :ref:`tutorial`.

In this tutorial, we will add further sectors to the electricity-only model from
:ref:`tutorial`, namely industry, transport, and buildings. This
requires processing of a few more raw data sources.

The sector-coupling code can be run as an overnight / greenfield scenario or
with multi-horizon investment with myopic foresight. Pathway analysis with
perfect foresight is under development. See also the documentation on
:ref:`foresight`.

Overnight Scenarios
===========================

Configuration
-------------

The default configuration file (``config/config.default.yaml``) is set up for running
overnight scenarios. Running a sector-coupled model unlocks many further
configuration options. In the example below, we say that the gas network should
be added and spatially resolved. We also say that the existing gas network may
be retrofitted to transport hydrogen instead.

.. literalinclude:: ../config/test/config.overnight.yaml
   :language: yaml
   :start-at: sector:
   :end-before: industry:

Documentation for all options will be added successively to :ref:`config`.

Scenarios can be defined like for electricity-only studies, but with additional configuration namespaces. Define scenario entries in ``config/scenarios.yaml`` and enable them via ``run.scenarios.enable: true`` to sweep different combinations of sector settings. See :doc:`wildcards` for the remaining dynamic placeholders.

Execution
---------

To run an overnight / greenfiled scenario with the specifications above, run

.. code:: console

    $ snakemake -call all --configfile config/test/config.overnight.yaml

Running this target orchestrates data retrieval, network preparation, and the final optimisation. Use ``snakemake -n all --configfile config/test/config.overnight.yaml`` to inspect the exact job summary for your setup; it will list the relevant ``build_*`` preprocessing tasks together with ``compose_network`` and ``solve_network``.

This covers the retrieval of additional raw data from online resources and
preprocessing data about the transport, industry, and heating sectors as well as
additional rules about geological storage and sequestration potentials, gas
infrastructure, and biomass potentials. The collection rule ``all`` will also
generate summary CSV files and plots after the network has been solved
successfully.



Inspect the generated DAG for your configuration with::

    snakemake --dag all --configfile config/test/config.overnight.yaml | dot -Tpng -o sector-dag.png

|

Myopic Foresight Scenarios
===================================

Configuration
-------------

To activate the myopic foresight mode, set

.. code:: yaml

    foresight: myopic

Scenarios can be defined like for electricity-only studies, but with additional
configuration namespaces. The myopic foresight mode relies on the top-level
``planning_horizons`` list to define the sequence of investment horizons:

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at: planning_horizons:
   :end-before: countries:

For allowed wildcard values (e.g., ``{horizon}``), refer to :ref:`wildcards`.

In the myopic foresight mode, you can tweak for instance exogenously given transition paths, like the one for
the share of primary steel production we change below:

.. literalinclude:: ../config/test/config.myopic.yaml
   :language: yaml
   :start-at: industry:
   :end-before: solving:

Documentation for all options will be added successively to :ref:`config`.

Execution
---------

To run a myopic foresight scenario with the specifications above, run

.. code:: console

    $ snakemake -call all --configfile config/test/config.myopic.yaml

This adds the sequential warm-start logic from :mod:`compose_network`, so the dry-run will list the same preprocessing tasks plus the additional per-horizon iterations. Use ``snakemake --dag all --configfile config/test/config.myopic.yaml`` to visualise the resulting graph.

|


Scaling-Up
==========

If you now feel confident and want to tackle runs with larger temporal, technological and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``all`` again without providing the test
configuration file.

.. code:: console

    $ snakemake -call purge
    $ snakemake -call all

.. note::

    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    .. code:: console

        $ snakemake -call all -n
