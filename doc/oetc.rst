.. SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
.. SPDX-FileCopyrightText: Open Energy Transition gGmbH
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
OETC Integration
##########################################

The OETC platform allows PyPSA-Eur to leverage cloud-based resources for
solving energy system models. This integration enables solving of optimization problems
using remote compute resources when local computational capacity is insufficient.

.. note::

   Using OETC is **optional** and not required to run PyPSA-Eur.

.. _oetc_overview:

Overview
========

OETC provides a cloud-based optimization service that can handle PyPSA networks by:

- Offloading optimization tasks to remote compute resources
- Supporting various solvers (Highs, Gurobi, etc.) on cloud infrastructure
- Automatically managing authentication and resource allocation
- Providing configurable compute resources (CPU cores, disk space)

The integration is implemented through the ``linopy.oetc`` module and is configured
within the PyPSA-Eur solving configuration.

.. _oetc_setup:

Setup and Authentication
========================

To use OETC with PyPSA-Eur, you need to:

1. **Set up environment variables** for authentication based on your OETC credentials:

   .. code:: bash

       export OETC_EMAIL="your-email@example.com"
       export OETC_PASSWORD="your-password"

2. **Configure OETC settings** in your configuration file (see :ref:`oetc_config`).

.. _oetc_config:

Configuration
=============

OETC integration is configured within the ``solving`` section of your configuration file.
Add an ``oetc`` subsection to enable cloud-based optimization:

.. code:: yaml

    solving:
      # ... other solving options ...
      oetc:
        name: "job-name" # Arbitrary human readable job identifier
        authentication_server_url: "https://auth.oetc.example.com"
        orchestrator_server_url: "https://orchestrator.oetc.example.com"
        compute_provider: "GCP"  # Currently only GCP is supported
        cpu_cores: 4 # This also sets the amount of RAM by a factor 8
        disk_space_gb: 20

.. note::

   The ``solver`` and ``solver_options`` are automatically passed from the main
   solving configuration and do not need to be specified in the OETC section.

.. _oetc_usage:

Usage
=====

Once configured, OETC integration works transparently with the standard PyPSA-Eur
workflow. If a oetc configuration is present, PyPSA-Eur will automatically use
the OETC platform for the model optimization during the solve network rule.

An example of a log entry when using OETC:
.. code::
    INFO:linopy.oetc:OETC - Signing in...
    INFO:linopy.oetc:OETC - Signed in
    INFO:linopy.oetc:OETC - Fetching user GCP credentials...
    INFO:linopy.oetc:OETC - Fetched user GCP credentials
    INFO:linopy.oetc:OETC - Submitting compute job...
    INFO:linopy.oetc:OETC - Compute job 992616b6-0f24-4ab5-8675-ba8098d795fa started
    INFO:linopy.oetc:OETC - Waiting for job 992616b6-0f24-4ab5-8675-ba8098d795fa to complete...
    INFO:linopy.oetc:OETC - Job 992616b6-0f24-4ab5-8675-ba8098d795fa status: PENDING, checking again in 30 seconds...
    INFO:linopy.oetc:OETC - Job 992616b6-0f24-4ab5-8675-ba8098d795fa status: PENDING, checking again in 45 seconds...
    INFO:linopy.oetc:OETC - Job 992616b6-0f24-4ab5-8675-ba8098d795fa completed successfully!
    INFO:linopy.oetc:OETC - Model solved successfully. Status: ok
    INFO:linopy.oetc:OETC - Objective value: 3.11e+07
