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

.. admonition:: What is the difference between snapshots and planning horizons?

.. admonition:: How can I resample an hourly time series to match my model's time resolution?

.. admonition:: In a myopic optimisation, how can I impose lower and upper bounds on renewable capacity expansion?

.. admonition:: How can I change the cost assumptions for natural gas, and can I use a time-varying gas price?

.. admonition:: Why does ``solve_network`` show the warning that the ``EU`` bus has no attached components?


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

.. admonition:: How can I model electricity imports from outside the model scope (e.g. an HVDC link to North Africa)?


Results & postprocessing
===============================================

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

