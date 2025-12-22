.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _faq_troubleshooting:

###############################################
FAQ and Troubleshooting
###############################################

This section contains answers to Frequently Asked Questions (FAQ) and common troubleshooting tips.

-----------------------------------------------

.. admonition:: FAQ
   :class: dropdown

   .. raw:: html

      <details>
      <summary><strong>I added x/y to the model, now it takes much longer to solve, why?</strong></summary>
      <p>
      Possible reasons include:
      <ul>
      <li>Adding new variables or constraints that narrow the feasible space significantly.</li>
      <li>Unfavourable bus clustering, leading to bottlenecks, e.g, very small clusters with high demand.</li>
      </p>
      </details>

-----------------------------------------------

.. admonition:: FAQ
   :class: dropdown

   .. raw:: html

      <details>
      <summary><strong>How can I build the documentation locally?</strong></summary>
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
      </details>