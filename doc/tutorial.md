<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Tutorial: Electricity-Only {#tutorial}

<iframe width="832" height="468" src="https://www.youtube.com/embed/mAwhQnNRIvs" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

!!! note
    If you have not done it yet, follow the [installation](installation.md) steps first.

In this tutorial, we will build a heavily simplified power system model for
Belgium. But before getting started with **PyPSA-Eur** it makes sense to be familiar
with its general modelling framework [PyPSA](https://pypsa.readthedocs.io).

Running the tutorial requires limited computational resources compared to the
full model, which allows the user to explore most of its functionalities on a
local machine. The tutorial will cover examples on how to configure and
customise the PyPSA-Eur model and run the ``snakemake`` workflow step by step
from network creation to the solved network. The configuration for the tutorial
is located at ``config/test/config.electricity.yaml``. It includes parts deviating from
the default config file ``config/config.default.yaml``. To run the tutorial with this
configuration, execute

```console
$ snakemake -call results/test-elec/networks/solved_2050.nc --configfile config/test/config.electricity.yaml
```

This configuration is set to download a reduced cutout via the rule `retrieve_cutout`.
For more information on the data dependencies of PyPSA-Eur, continue reading [Retrieving Data](retrieve.md#data).

## How to configure runs?

The model can be adapted to only include selected countries (e.g. Belgium) instead of all European countries to limit the spatial scope.

```yaml
{{ yaml_section("countries", source="test/config.electricity.yaml") }}
```

Likewise, the example's temporal scope can be restricted (e.g. to a single week).

```yaml
{{ yaml_section("snapshots", source="test/config.electricity.yaml") }}
```

It is also possible to allow less or more carbon-dioxide emissions. Here, we limit the emissions of Belgium to 100 Mt per year.

```yaml
{{ yaml_section("co2_budget", source="test/config.electricity.yaml") }}
```

PyPSA-Eur also includes a database of existing conventional powerplants.
We can select which types of existing powerplants we like to be extendable:

```yaml
{{ yaml_section("electricity.extendable_carriers", source="test/config.electricity.yaml") }}
```

To accurately model the temporal and spatial availability of renewables such as
wind and solar energy, we rely on historical weather data. It is advisable to
adapt the required range of coordinates to the selection of countries.

```yaml
{{ yaml_section("atlite", source="test/config.electricity.yaml") }}
```

We can also decide which weather data source should be used to calculate
potentials and capacity factor time-series for each carrier. For example, we may
want to use the ERA-5 dataset for solar and not the default SARAH-3 dataset.

```yaml
{{ yaml_section("renewable.solar.cutout") }}
```

Finally, it is possible to pick a solver. For instance, this tutorial uses the
open-source solver HiGHS.

```yaml
{{ yaml_section("solving.solver", source="test/config.electricity.yaml") }}
```

Note, that ``config/test/config.electricity.yaml`` only includes changes relative to
the default configuration. There are many more configuration options, which are
documented at [Configuration](configuration.md).

### Directory Structure and Configuration Settings

It's important to understand how certain configuration settings affect the directory structure in PyPSA-Eur:

- ``run.name`` determines the subdirectory within the ``results`` folder (e.g., ``results/test-elec/networks/...``)
- ``run.shared_resources.policy`` determines the subdirectory within the ``resources`` folder (e.g., ``resources/test/networks/...``)

These settings work together to organize model runs:

- Final model outputs are always stored in ``results/[run.name]/...``
- Intermediate files can be either:
    - Specific to a run: ``resources/[run.shared_resources.policy]/...`` (if policy is a string)
    - Shared between runs: ``resources/...`` (if policy is ``true``)
    - Not shared between runs: ``resources/[run.name]`` (if policy is ``false``)
    - Partially shared: If policy is ``"base"``, some common files are shared while others remain run-specific

For this tutorial, with ``run.name: "test-elec"`` and ``run.shared_resources.policy: "test"``,
intermediate resources are stored in ``resources/test/...`` while results are in ``results/test-elec/...``.

The implementation of this behavior can be found in ``scripts/_helpers.py``.

## How to use ``snakemake`` rules?

Open a terminal, go into the PyPSA-Eur directory, and activate the ``pypsa-eur`` environment with

```console
$ pixi shell
```

Let's say based on the modifications above we would like to solve a very simplified model
clustered down to 6 buses and every 24 hours aggregated to one snapshot. The command

```console
$ snakemake -call results/test-elec/networks/solved_2050.nc --configfile config/test/config.electricity.yaml
```

orders ``snakemake`` to run the rule [solve_network][] that produces the solved network and stores it in ``results/test-elec/networks`` with the name ``solved_2050.nc``:

```python
rule solve_network:
    input:
        network=resources("networks/composed_{horizon}.nc"),
    output:
        network=RESULTS + "networks/solved_{horizon}.nc",
    log:
        solver=normpath(RESULTS + "logs/solve_network/solver_{horizon}.log"),
        memory=RESULTS + "logs/solve_network/memory_{horizon}.log",
        python=RESULTS + "logs/solve_network/python_{horizon}.log",
    benchmark:
        (RESULTS + "benchmarks/solve_network_{horizon}.log")
    shadow:
        shadow_config
    threads: solver_threads
    resources:
        mem_mb=config_provider("solving", "mem_mb"),
        runtime=config_provider("solving", "runtime", default="6h"),
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        planning_horizons=config_provider("planning_horizons"),
        sector=config_provider("sector"),
        co2_sequestration_potential=config_provider(
            "sector", "co2_sequestration_potential"
        ),
        custom_extra_functionality=input_custom_extra_functionality,
    script:
        "../scripts/solve_network.py"
```

This triggers a workflow of multiple preceding jobs that depend on each rule's inputs and outputs.
Inspect the full directed acyclic graph with:

```console
$ pixi run dot -c && snakemake --dag results/test-elec/networks/solved_2050.nc --configfile config/test/config.electricity.yaml | dot -Tpng -o dag.png
```

In the terminal, this will show up as a list of jobs to be run:

```console
Building DAG of jobs...
Job stats:
job                                      count
-------------------------------------  -------
add_transmission_projects_and_dlr            1
base_network                                 1
build_electricity_demand                     1
build_electricity_demand_base                1
build_line_rating                            1
build_monthly_prices                         1
build_osm_boundaries                         4
build_population_layouts                     1
build_powerplants                            1
build_renewable_profiles                     3
build_shapes                                 1
build_solar_rooftop_potentials               1
build_solar_thermal_profiles                 1
build_transmission_projects                  1
cluster_network                              1
compose_network                              1
determine_availability_matrix                3
process_cost_data                            1
retrieve_cost_data                           1
retrieve_cutout                              1
retrieve_databundle                          1
retrieve_eez                                 1
retrieve_electricity_demand                  1
retrieve_eurostat_data                       1
retrieve_jrc_ardeco                          1
retrieve_monthly_co2_prices                  1
retrieve_monthly_fuel_prices                 1
retrieve_nuts_2021_shapes                    1
retrieve_osm_boundaries                      4
retrieve_osm_prebuilt                        1
retrieve_synthetic_electricity_demand        1
retrieve_worldbank_urban_population          1
simplify_network                             1
solve_network                                1
time_aggregation                             1
total                                       46
```


``snakemake`` then runs these jobs in the correct order.

A job (here ``build_powerplants``) will display its attributes and normally some logs below this block:

```console
rule build_powerplants:
    input: resources/test-elec/networks/clustered.nc, data/custom_powerplants.csv
    output: resources/test-elec/powerplants.csv
    log: logs/test-elec/build_powerplants.log
    jobid: 40
    benchmark: benchmarks/test-elec/build_powerplants
    reason: Missing output files: resources/test-elec/powerplants.csv
    resources: tmpdir=<TBD>, mem_mb=7000, mem_mib=6676
```

Once the whole worktree is finished, it should state so in the terminal.

You will notice that many intermediate stages are saved, namely the outputs of each individual ``snakemake`` rule.

You can produce any output file occurring in the ``Snakefile`` by running

```console
$ snakemake -call <output file>
```

For example, you can explore the evolution of the PyPSA networks by running

1. ``snakemake -call resources/test-elec/networks/base.nc --configfile config/test/config.electricity.yaml``
2. ``snakemake -call resources/test-elec/networks/simplified.nc --configfile config/test/config.electricity.yaml``
3. ``snakemake -call resources/test-elec/networks/clustered.nc --configfile config/test/config.electricity.yaml``
4. ``snakemake -call resources/test-elec/networks/composed_2050.nc --configfile config/test/config.electricity.yaml``
5. ``snakemake -call results/test-elec/networks/solved_2050.nc --configfile config/test/config.electricity.yaml``

To run all scenario combinations defined in ``config/scenarios.yaml`` (enable this via ``run.scenarios.enable: true``),
you can use the collection rule ``solve_networks``.

```console
$ snakemake -call solve_networks --configfile config/test/config.electricity.yaml
```

If you now feel confident and want to tackle runs with larger temporal and
spatial scope, clean-up the repository and after modifying the ``config/config.yaml`` file
target the collection rule ``solve_networks`` again without providing the test
configuration file.

```console
$ snakemake -call purge
snakemake -call solve_networks
```

!!! note
    It is good practice to perform a dry-run using the option `-n`, before you
    commit to a run:

    ```console
    $ snakemake -call solve_networks -n
    ```

## How to analyse results?

The solved networks can be analysed just like any other PyPSA network (e.g. in
Jupyter Notebooks).

```python
import pypsa

n = pypsa.Network("results/test-elec/networks/solved_2050.nc")
```

For inspiration, read the [examples section in the PyPSA documentation](https://pypsa.readthedocs.io/en/latest/examples-basic.html).
