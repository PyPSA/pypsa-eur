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

- `Installation of PuTTy `<https://www.ibm.com/products/ilog-cplex-optimization-studio>`_ ## not sure about that one


Step 1 - Google Cloud Platform registration
====================

First of all, you have to register at the `Google Cloud Platform <https://console.cloud.google.com>`_ (GCP). No magic.
You only need an active bank account. Don't worry they won't charge you for anything as far as you are in the free trial budget and don't want to continue the GCP.
( Top secret tip. If you run out of your first 300$ free trial budget, you can simply ask a friend or family member to set up another account for you. Though, we don't recommend that as long-term solution.)


Step 2 - Create your Virtual Machine instance
===============================

With the following steps we create a Virtual Machine (VM) on Google Cloud.

- Click on the `GCP Dashboard <https://console.cloud.google.com/home/dashboard>`_.
- Click at the "COMPUTE" header, on the "Compute Engine" and then on the "VM instance".
- Click on create.
- Click on new VM instance.

Now the a window with the machine detail will open. You have to configure the following things:

- Name. Set a name for your VM (Chose your name well i.e. "climate-casino". You cannot edit the name after saving the instance.)
- Region. You can keep us-central1 (Iowa), since it is a 'cheap' computational region. Sometimes your machine is limited in a specific region, dont worry to pick another region.
- Machine configuration. The machine configuration sets how powerful your VM is. For the set-up stage we suggest that you using a 1 vCPU and 3.75 GB memory, N1 series machine. We suggest that because every operating second cost money (Your valuable free trial money!). In a later stage you can edit your machine configuration without problems. So use a cheap machine type configuration to handle/ transfer data and only when everything is ready and tested, your expensive machine type, for instance a 8 vCPU with 160 GB memory (check "snakemake -j -n 1 solve_all_elec_networks" as a dry run to see what PyPSA-Eur or PyPSA-Eur-Sec requires for memory. Usually PyPSA can't handle not more than 8 vCPU so no more cores than 8 are required. But the memory requirements can often vary depending on the spatial and time detail of your simulation (I.e. we computed an hourly, 181 node full EU-network, with 8 vCPU and 150 GB memory since the dry-run gave as 135 GB memory as min. requirement.)
- Boot disk. As default, your VM is created with 10 GB. Depending on how much you want to handle on one VM you should increase the disk size. We recommend a disk size of 100 GB for a safe start (cost roughly 8$ per month), the disk can be resized at any later stage as additional disk.

- Click on create and celebrate your first VM on GCP.

Installation of Cloud SDK
=========================

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
