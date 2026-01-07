.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _faq_troubleshooting:

###############################################
FAQ and Troubleshooting
###############################################

This section contains answers to Frequently Asked Questions (FAQ) and common troubleshooting tips.

-----------------------------------------------


General
===============================================

.. admonition:: I'm having trouble installing PyPSA-Eur or getting started. Where should I start?

   The most common installation issues involve Python environment setup and solver configuration. We recommend using ``pixi`` for environment management. For solver setup, HiGHS is included by default for testing, but commercial solvers are supported as well. See :doc:`installation` for detailed platform-specific instructions, solver configuration guidance and alternative environment manager if you prefer using ``conda`` (legacy).

.. admonition:: My workflow is failing or producing unexpected results. How do I troubleshoot?

   Start by running ``snakemake -call -n`` (dry-run) to validate workflow structure without execution. Then, check log files in ``logs/`` and verify intermediate results at each workflow stage. For persistent issues, see :doc:`support` for community assistance channels.

.. admonition:: I would like to develop and maintain a custom fork of PyPSA-Eur with additional features. How can I best manage this?

   We recommend forking the PyPSA-Eur repository on GitHub and creating a dedicated branch for your custom features. Regularly sync your fork with the upstream repository to incorporate updates and improvements. Use feature branches for individual changes and consider submitting pull requests to contribute back to the main project if applicable.

.. admonition:: What is the difference between snapshots and planning horizons?

   **Snapshots** define the temporal resolution of system operation. They represent the individual timesteps (e.g. hourly, 4-hourly, or segments of individual lengths) over which the dispatch of technologies is optimised.

   **Planning horizons**, used for example in overnight or myopic optimisation workflows, define distinct investment periods. They either represent an individual year (e.g. 2030) or divide the long-term investment period into multiple stages (e.g. 2030, 2040, and 2050). For myopic optimisations, investments made in one horizon influence the next.

   In short: snapshots handle *operational time*, planning horizons handle *investment time*.

.. admonition:: In a myopic optimisation, how can I impose lower and upper bounds on renewable capacity expansion?

.. admonition:: How can I change the cost assumptions for natural gas, and can I use a time-varying gas price?


Solving
===============================================

.. admonition:: My model takes much longer to solve after adding or modifying X. Why?

   A noticeable increase in solve time usually means that the optimisation problem has become more complex. This may occur when additional variables or constraints restrict the feasible space, when the bus clustering introduces bottlenecks (for example, very small regions with high demand), or when the temporal or spatial resolution is set too high.

   For reference, PyPSA-Eur typically solves efficiently with 128 spatial clusters and 4-hour temporal resolution, assuming no custom constraints are added. If you increase either the temporal or spatial resolution, you may need to decrease the other to keep solution times manageable.

.. admonition:: What can I do if my dispatch model reports “infeasible or unbounded”?

   This message typically indicates that the model cannot satisfy demand in all
   timesteps or that some generators have lower bounds exceeding the available
   load. If not already enabled, activate ``load_shedding`` and
   `curtailment mode <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#solving>`__ in the solving configuration. These options introduce a high-cost load-shedding generator and a curtailment
   generator, allowing the optimisation to remain feasible.

   If enabling these options resolves the issue, examine the resulting dispatch to
   determine the bus regions and snapshots where the system was unable to meet
   demand. A simple way to quantify the total load shedding is:

   .. code-block:: python

      load_shedding_i = n.generators[n.generators.carrier == "load"].index
      load_shedding = (
          n.generators_t.p[load_shedding_i]
          .mul(n.snapshot_weightings.generators, axis=0)
          .sum()
      )
      print(load_shedding)

.. admonition:: Why does ``solve_network`` show the warning that the ``EU`` bus has no attached components?

   The ``EU`` bus is introduced as an aggregate carrier bus used to track system-wide quantities—such as total fuel use (e.g. hydrogen, biomass, CO₂) across all regions. This is particularly useful when treating the respective fuel in a copperplated configuration. Because this bus is not intended to represent a physical location and no components are directly attached to it in the optimisation problem, PyPSA
   emits a warning. The warning is harmless and can safely be ignored.


Model configuration & customisation
===============================================

.. admonition:: How can I change the spatial resolution of my model?

   Several configuration options control the spatial resolution in PyPSA-Eur. Available mechanisms:

   - Adjust the number of regions via the `clusters <https://pypsa-eur.readthedocs.io/en/latest/wildcards.html#clusters>`__ wildcard.
   - Restrict the geographical scope using the `countries <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#countries>`__ configuration.
     Note that the number of clusters must exceed the number of countries.
   - Configure how regions are created through the `clustering <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#clustering>`__ section,
     including NUTS-based clustering (``administrative``) or custom weightings
     (``focus_weights``).
   - Provide a fully custom spatial clustering by setting ``mode: custom_busshapes`` and
     placing your geometry file at
     ``data/busshapes/base_s_{w.clusters}_{base_network}.geojson``.

.. admonition:: Why do some countries (e.g. GB or DK) have at least two nodes even when I specify only one cluster per country?

   Some countries span multiple synchronous zones. Examples include Denmark (DK1 and DK2) and Great Britain versus Northern Ireland. To represent these electrical boundaries correctly, PyPSA-Eur assigns one bus region per synchronous zone, even if the configuration requests only a single cluster for the whole country.

   If you prefer a purely administrative grouping without splitting along synchronous zone boundaries, set the ``clustering`` option to use the `NUTS classification <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#clustering>`__ (``administrative``) in the configuration.

.. admonition:: How can I ensure that my custom changes are correctly implemented in the model?

   For debugging and tracing, it is a good practice to implement logs (``logging`` module in Python) in your custom scripts. This allows you to track the flow of execution and identify where things might be going wrong.

   Possible steps: Use ``n.statistics.energy_balance`` (`docs <https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.energy_balance>`__) to ensure that the technology you added is expanded and dispatched. Use ``n.statistics.opex`` (`docs <https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.opex>`__) and ``n.statistics.capex`` (`docs <https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.capex>`__) to ensure that the overall system costs stay in a reasonable range compared to the previous (or counterfactual) run.

.. admonition:: How can I resample an hourly time series to match my model's time resolution?

   Use PyPSA's built-in `resampling functions <https://docs.pypsa.org/latest/api/other/common/#pypsa.common.resample_timeseries>`__, which correctly aggregates time series and preserves all required metadata. For example, to convert an hourly profile to 4-hour resolution:

   .. code-block:: python

      from pypsa.common import resample_timeseries

      ts_resampled = resample_timeseries(df, freq="4H")

   Please be aware that the dataframe's index must be a DateTimeIndex.

.. admonition:: How can I model electricity imports from outside the model scope (e.g. an HVDC link to North Africa)?


Results & postprocessing
===============================================

.. admonition:: How can I start analysing my results (e.g. energy balances) after a successful model run?

   PyPSA-Eur already generates a rich set of default outputs. These are stored in the ``results/`` directory and include maps, plots, and summary tables that provide a first overview of capacities, dispatch patterns, prices, and line loadings. Reviewing these default figures is often the quickest way to understand the high-level behaviour of your scenario.

   For more detailed analysis, use PyPSA's `built-in statistics API <https://docs.pypsa.org/latest/api/networks/statistics/>`__, which offers a variety of functions to extract and summarise key performance indicators. These statistics can be further processed using standard Python libraries such as Pandas.

   .. code-block:: python

      eb = n.statistics.energy_balance()
      print(eb)

   You can filter statistics by carrier, component, bus, or region using PyPSA's flexible filtering API. See the `documentation <https://docs.pypsa.org/stable/user-guide/statistics/#filtering>`__ for more details:

   Other useful entry points include:

   .. code-block:: python

      n.statistics.capex()      # investment costs
      n.statistics.opex()       # operational costs
      n.statistics.generator()  # generator-level summaries

   These tools help you go beyond the default outputs and systematically explore the solved system.

.. admonition:: I observe a very high marginal electricity price in a single timestep. Is this normal?

   This can be normal, but it is not always expected. In capacity expansion models with perfectly inelastic demand (the PyPSA-Eur default), a very high marginal price can appear when the system is operating at its limits. In such a timestep, supplying an additional unit of demand would require violating a constraint (for example, insufficient generation, storage, or transmission) or would trigger load shedding. The corresponding shadow price of the demand balance constraint then becomes large.

   However, unusually high prices can also signal modelling issues such as missing flexibility options, too few clusters causing artificial bottlenecks, an overly restrictive constraint, or incorrectly configured capacity limits. For an illustration of why marginal prices can spike when demand is perfectly inelastic, see the `demand-elasticity example <https://docs.pypsa.org/latest/examples/demand-elasticity/#perfectly-inelastic-demand-up-to-voll>`__ in the PyPSA documentation.


Contribution & support
===============================================

.. admonition:: How do I compile the documentation locally?
   :class: tip

   The documentation is built with `Sphinx <https://www.sphinx-doc.org/>`__
   using the dedicated ``doc`` Pixi environment. You can build it either through
   the predefined Pixi task or manually with Sphinx.

   **Option 1 — Use the Pixi task (recommended)**

   .. code-block:: bash

      pixi run -e doc build-docs _build html

   This creates the HTML output in ``doc/_build/html``.  
   Open ``doc/_build/html/index.html`` in your browser to view the result.

   **Option 2 — Build manually with Sphinx**

   .. code-block:: bash

      pixi install -e doc
      pixi shell -e doc
      cd doc
      make html

   This also generates the HTML documentation in ``doc/_build/html``.

