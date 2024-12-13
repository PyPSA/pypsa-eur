<!--
SPDX-FileCopyrightText: 2017-2024 The PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/pypsa/pypsa-eur?include_prereleases)
[![Test workflows](https://github.com/pypsa/pypsa-eur/actions/workflows/test.yaml/badge.svg)](https://github.com/pypsa/pypsa-eur/actions/workflows/test.yaml)
[![Documentation](https://readthedocs.org/projects/pypsa-eur/badge/?version=latest)](https://pypsa-eur.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/pypsa/pypsa-eur)
[![Zenodo PyPSA-Eur](https://zenodo.org/badge/DOI/10.5281/zenodo.3520874.svg)](https://doi.org/10.5281/zenodo.3520874)
[![Zenodo PyPSA-Eur-Sec](https://zenodo.org/badge/DOI/10.5281/zenodo.3938042.svg)](https://doi.org/10.5281/zenodo.3938042)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.14.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![REUSE status](https://api.reuse.software/badge/github.com/pypsa/pypsa-eur)](https://api.reuse.software/info/github.com/pypsa/pypsa-eur)
[![Stack Exchange questions](https://img.shields.io/stackexchange/stackoverflow/t/pypsa)](https://stackoverflow.com/questions/tagged/pypsa)

# EU-Flex Platform
<img src="https://raw.githubusercontent.com/open-energy-transition/oet-website/main/assets/img/oet-logo-red-n-subtitle.png" alt="Open Energy Transition Logo" width="260" height="100" align="right">


## Overview
The EU-Flex Platform is an implementation of ACER's EU-wide flexibility assessment platform based on PyPSA-Eur. This project develops a comprehensive modeling tool for assessing power system flexibility needs across the European Union, as mandated by Regulation (EU) 2024/1747.

## Project Scope
This repository covers Phase 1 (Implementation) and Phase 2 (Testing) of the EU-Flex Platform development:

### Phase 1: Implementation (Dec 2024 - Nov 2025)
- Development of core modeling capabilities using PyPSA
- Implementation of flexibility assessment functionalities
- Integration of various data sources and power system components
- Development of APIs and user interfaces
- Setting up security and operational protocols

### Phase 2: Testing (Dec 2025 - June 2026)
- Comprehensive testing with complete datasets
- Performance optimization and refinement
- Integration testing with real-world scenarios
- System validation and quality assurance
- Methodology alignment and adjustments

## Key Features (WIP)
- EU-wide power system modeling
- Flexibility needs assessment across multiple timeframes
- Integration of renewable energy sources
- Analysis of demand-side response and storage options
- Cross-border capacity evaluation
- Stochastic analysis capabilities

## Repository structure

* `benchmarks`: will store `snakemake` benchmarks (does not exist initially)
* `config`: configurations used in the study
* `cutouts`: will store raw weather data cutouts from `atlite` (does not exist initially)
* `data`: includes input data that is not produced by any `snakemake` rule
* `doc`: includes all files necessary to build the `readthedocs` documentation of PyPSA-Eur
* `envs`: includes all the `mamba` environment specifications to run the workflow
* `logs`: will store log files (does not exist initially)
* `notebooks`: includes all the `notebooks` used for ad-hoc analysis
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `rules`: includes all the `snakemake`rules loaded in the `Snakefile`
* `resources`: will store intermediate results of the workflow which can be picked up again by subsequent rules (does not exist initially)
* `results`: will store the solved PyPSA network data, summary files and plots (does not exist initially)
* `scripts`: includes all the Python scripts executed by the `snakemake` rules to build the model

## Technology Stack
- [OET/PyPSA-Eur](https://github.com/open-energy-transition/pypsa-eur) used via the [OET's soft-fork strategy](https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy). The model is described in the [documentation](https://pypsa-eur.readthedocs.io) and in the paper
[PyPSA-Eur: An Open Optimisation Model of the European Transmission System](https://arxiv.org/abs/1806.01613), 2018, [arXiv:1806.01613](https://arxiv.org/abs/1806.01613).
- [Snakemake](https://snakemake.readthedocs.io) for workflow management

## Installation and Usage

### 1. Installation

Clone the repository:

    git clone https://github.com/open-energy-transition/{{repository}}

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Users may also prefer to use [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) or [conda](https://docs.conda.io/projects/conda/en/latest/index.html). Using `mamba`, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

Activate the newly created `{{project_short_name}}` environment:

    mamba activate {{project_short_name}}

### 2. Run the analysis

    snakemake -call

This will run all analysis steps to reproduce results and build the report.

To generate a PDF of the dependency graph of all steps `resources/dag.pdf` run:

    snakemake -c1 dag

<sup>*</sup> Open Energy Transition (g)GmbH, Königsallee 52, 95448 Bayreuth, Germany

## Access and Usage
This platform is being developed for ACER (European Union Agency for the Cooperation of Energy Regulators). Access and usage rights are restricted according to ACER's policies and regulations.
