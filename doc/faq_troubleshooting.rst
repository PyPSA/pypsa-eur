.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _faq_troubleshooting:

###############################################
FAQ and Troubleshooting
###############################################

This section contains answers to Frequently Asked Questions (FAQ) and common troubleshooting tips.

-----------------------------------------------

.. raw:: html

   <h4>❓I added x/y to the model, now it takes much longer to solve, why?</h4>
      Possible reasons include:
      <li>Adding new variables or constraints that narrow the feasible space significantly.</li>
      <li>Unfavourable bus clustering, leading to bottlenecks, e.g, very small clusters with high demand.</li>
      <li>High temporal and/or spatial resolution.</li>
      </p>
      <p>
      How trouble shoot:
      <li>Without custom constraints and variables pypsa-eur usually solves with a spatial resolution of 128 clusters and a time resolution of 4H.</li>
      <li>If you increase either the temporal or spatial resolution you might have to decrease the other one.</li>
      </p>

   <h4>❓ I changed a parameter or added a technology. How can I check that I did it correctly?</h4>
      <li>Ideally you compare your new results with a standard pypsa-eur model.
      <li>Use <a href="https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.energy_balance" ><code>n.statistics.energy_balance</code></a> to ensure that the technology you added is expanded and dispatched.
      <li>Use <a href="https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.opex"><code>n.statistics.opex</code></a> + <a href="https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.capex"><code>n.statistics.capex</code></a> to ensure that the overall system costs stay in a reasonable range compared to the standard run.

   <h4>❓ I want to change the spatial resolution of my model. What options do I have?</h4>
      <li>The config file offers a number of possibilities how to manipulate the spatial resolution
      <li>The wildcard <a href="https://pypsa-eur.readthedocs.io/en/latest/wildcards.html#clusters"><code>clusters</code></a> allows a quick adjustment of the number of regions
      <li>The config <a href="https://pypsa-eur.readthedocs.io/en/latest/configuration.html#countries"><code>countries</code></a> specifies which countries are included in the model. Note that the number of clusters must be greater than the number of countries.
      <li>The config <a href="https://pypsa-eur.readthedocs.io/en/latest/configuration.html#clustering"><code>clustering</code></a> gives a number of possibilities on how to distribute the regions. There is the option to use NUTS classification <code>administrative</code> or assign weights for each of the countries <code>focus_weights</code>.
      <li>You can specify your own clustering by setting <code>mode: custom_busshapes</code> and store the corresponding file under <code>data/busshapes/base_s_{w.clusters}_{base_network}.geojson</code>.
   
   <h4>❓ I'm running a dispatch model and get the status "infeasible or unbounded". What can I do?</h4>
      <li>Most likely the load cannot be met in all timesteps or there is a lower bound for the generation which is higher than the demand.
      <li>If not activated, enable <a href="https://pypsa-eur.readthedocs.io/en/latest/configuration.html#solving"><code>load_shedding</code></a> and <code>curtailment_mode</code>.
      <li>This adds a generator with high costs for load shedding and a generator that allows to curtail electricity.
      <li>If this solved your numerical trouble, you can investigate the dispatch of the generators to isolate the bus region and snapshot in which the system cannot meet the load.
      </p>
      <pre>
      load_shedding_i = n.generators[n.generators.carrier=="load"].index
      load_shedding = n.generators_t.p[load_shedding_i].mul(n.snapshot_weightings.generators, axis=0).sum()
      print(load_shedding)
      </pre>
      <p>

   <h4>❓ When I analyze my model I notice that the marginal electricity price is very high in a single timestep. Is that normal?</h4>
      <li>This is normal behavior of an energy system model.
      <li>There is a timestep in which an additional load increment would require capacity expansion which leads to the very high marginal price.
      <li>See also the <a href="https://docs.pypsa.org/latest/examples/demand-elasticity/#perfectly-inelastic-demand-up-to-voll">example for demand elasticity</a>

   <h4>❓ Why are there at least two nodes for certain countries (GB/DK) even though I specified only one node per country?</h4>
      <li>There are countries with different synchronous zones like Denmark and Great Britain.
      <li>To account for those, PyPSA-Eur is allocating one bus region per synchronous zone.
      <li>You can avoid that by setting <a href="https://pypsa-eur.readthedocs.io/en/latest/configuration.html#clustering"><code>clustering</code></a> to match the NUTS classification <code>administrative</code>

   <h4>❓ What is the difference between snapshots and planning horizons?</h4>

   <h4>❓ I have a hourly timeseries for a component that I want to resample to match my time resolution. Is there an easy way to do that?</h4>

   <h4>❓ How do I add electricity import across the model's scope e.g. via a HVDC link to North Africa?</h4>

   <h4>❓ In my myopic optimization, how can I add a lower and upper corridor for the build out of renewable energies?</h4>

   <h4>❓ I want to change the assumptions for the costs of natural gas. How can I do that and can I even pass a temporally changing price for gas?</h4>

   <h4>❓ The rule solve_network returns the warning "WARNING:pypsa.consistency:The following buses have no attached components, which can break the lopf: {'EU'}". What can I do about that?
      <li>The EU bus is added for plotting reasons which is why you can ignore that warning.

   <h4>❓ How can I build the documentation locally?</h4>
      <p>
      Activate the Pixi documentation environment and run Sphinx:
      </p>
      <pre>
      pixi install -e doc
      pixi shell -e doc
      cd doc
      make html
      </pre>
      <p>
      Then open <code>doc/_build/html/index.html</code>.
      </p>