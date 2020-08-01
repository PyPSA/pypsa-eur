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

Step 1 - Google Cloud Platform registration
====================

First of all, you have to register at the `Google Cloud Platform <https://console.cloud.google.com>`_ (GCP). No magic.
You only need an active bank account. Don't worry they won't charge you for anything as far as you are in the free trial budget and don't want to continue the GCP.
(Tip. If you run out of your first 300$ free trial budget, you can simply ask a friend or family member to set up another account for you. Though, we don't recommend that as long-term solution.)


Step 2 - Create your Virtual Machine instance
===============================

With the following steps we create a Virtual Machine (VM) on Google Cloud.

- Click on the `GCP Dashboard <https://console.cloud.google.com/home/dashboard>`_.
- Click at the "COMPUTE" header, on the "Compute Engine" and then on the "VM instance".
- Click on create.
- Click on new VM instance.

Now a window with the machine details will open. You have to configure the following things:

- Name. Set a name for your VM (Chose your name well i.e. "climate-casino". You cannot edit the name after saving the instance.)
- Region. You can keep us-central1 (Iowa), since it is a 'cheap' computational region. Sometimes your machine is limited in a specific region, dont worry to pick another region.
- Machine configuration. The machine configuration sets how powerful your VM is. For the set-up stage we suggest that you using a 1 vCPU and 3.75 GB memory, N1 series machine. We suggest that because every operating second cost money (Your valuable free trial money!). In a later stage you can edit your machine configuration without problems. So use a cheap machine type configuration to handle/ transfer data and only when everything is ready and tested, your expensive machine type, for instance a 8 vCPU with 160 GB memory (check "snakemake -j -n 1 solve_all_elec_networks" as a dry run to see what PyPSA-Eur or PyPSA-Eur-Sec requires for memory. Usually PyPSA can't handle not more than 8 vCPU so no more cores than 8 are required. But the memory requirements can often vary depending on the spatial and time detail of your simulation (I.e. we computed an hourly, 181 node full EU-network, with 8 vCPU and 150 GB memory since the dry-run gave as 135 GB memory as min. requirement.)
- Boot disk. As default, your VM is created with 10 GB. Depending on how much you want to handle on one VM you should increase the disk size. We recommend a disk size of 100 GB for a safe start (cost roughly 8$ per month), the disk can be resized at any later stage as additional disk.

- Click on create and celebrate your first VM on GCP.

Step 3 - Installation of Cloud SDK
===================================

- Download Google Cloud SDK `SDK <https://cloud.google.com/sdk>`_. Check that you are logged in in your google account. The link should lead you to the windows installation of Google Cloud SDK.
- Follow the "Quickstart for Windows - Before you begin" steps.
- After the successfull installation, close the Google Cloud SDK reopen it again. Type the following command into the "Google Cloud SDK Shell":

    .. code:: bash
        
        gcloud compute ssh [Your VM instance name] -- -L 8888:localhost:8888
        
- This command above will open a PuTTy command window that is connected to your Virtual Machine. Time to celebrate if it works!
- Now install all necessary tools. As little help, the first steps: 

    .. code:: bash
        
        sudo apt-get update
        sudo apt-get install bzip2 libxml2-dev
        sudo apt-get install wget
        wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh    (Check the http link. To be up to date with anaconda, check the Anaconda website https://www.anaconda.com/products/individual )
        ls  (to see what anaconda file to bash)
        bash Anaconda3-2020.07-Linux-x86_64.sh  
        source ~/.bashrc  
        
- Close and reopen the PuTTy file (-> open Google Cloud SDK -> initialize again with the command above to open the PuTTY command window). Now Conda can be listed with 'conda list'. Follow now the standard PyPSA/PyPSA-Eur/PyPSA-Eur-Sec installation pages to finalize your machine for any PyPSA modelling tasks (don't forget the solvers - for bigger simulations use commercial solvers such as Gurobi).
        
Step 4 - Installation of WinSCP
===================================  

For smooth data exchange between the VM and your local computer we recommend WinSCP.
Make sure that your instance is operating for the next steps.

- Download `WinSCP <https://winscp.net/eng/download.php>`_ and follow the default installation steps.
- Open WinSCP after the installation. A login window will open.
- Keep SFTP as file protocol.
- As host name insert the External IP of your VM (click in your internet browser on your GCP VM instance to see the external IP) 
- Set the User name in WinSCP to the name you see in your PuTTy window (check step 3 - for instance [username]@[VM-name]:~$)
- Click on the advanced setting. SSH -> Authentication. 
- Option 1. Click on the Tools button and "Install Public Key into Server..". Somewhere in your folder structure must be a public key. I found it with the following folder syntax on my local windows computer -> :\Users\Max\.ssh (there should be a PKK file). 
- Option 2 (alternative). Click on the Tools button and "Generate new key pair...". Save the private key at a folder you remember and add it to the "private key file" field in WinSCP. Upload the public key to the metadeta of your instance. 
- Click ok and save. Then click Login. If successfull WinSCP will open on the left side your local computer folder structure and on the right side the folder strucutre of your VM. (If you followed Option 2 and its not initially working. Stop your instance, refresh the website, reopen the WinSCP field. Afterwards your your Login should be successfull)

If you had struggle with the above steps, you could also try `this video <https://www.youtube.com/watch?v=lYx1oQkEF0E>`_.

.. note::
    Double check the External IP of your VM before you try to login with WinSCP. It's often a cause for an error.
..


Step 5 - Extra. Copying your instance with all its data/ paths included.
========================================================================
Especially if you think about operating several instance for quicker simulations, you can create a so called `"image" <https://console.cloud.google.com/compute/images?authuser=1&project=exalted-country-284917>`_ of the virtual machine. The "image" include all the data and software set-ups from your VM. Afterwards you can create a VM from an image and avoid all the installation steps above. 

Important points when to solve networks.
========================================================================
If you use the GCP with the default PyPSA-Eur settings, your free budget will dissapear quickly. The following tips should help you to be efficient with your free budget.

- Test always in low resolution networks. I.e single country with 5 nodes, 24h time resolution for 2 month data.
- Adjust your solver in the config.yaml file. 
- 1. At "solving:" reset "skip_iterations:" from "false" to "true". This will lead to a single solver iteration which is often precise enough, since the following iteration barely change the objective value. 
- 2. At "solver:" and "FeasibilityTol:" increase the tolerance to 1.e-4. This will lead to a slightly less accurate objective value but lowers drastically the computational requirements.
