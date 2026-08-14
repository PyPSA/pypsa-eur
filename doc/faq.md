<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# FAQ and Troubleshooting

This section contains answers to Frequently Asked Questions (FAQ) and common troubleshooting tips.

## General

??? note "Snakemake `retrieve_*` rules fail with `SSL: CERTIFICATE_VERIFY_FAILED` error"
    If your workflow fails while downloading input data (for example during `snakemake -n`), and you may see errors like this:

    ```console
    $ pixi run snakemake -n --cores 1
    Building DAG of jobs...
    ConnectError while get'ing https://data.pypsa.org/workflows/eur/corine/v18_5/manifest.yaml
    [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self-signed certificate in certificate chain (_ssl.c:1032)
    ...
    httpx.ConnectError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self-signed certificate in certificate chain (_ssl.c:1032)
    ```

    This typically happens on networks where a university/company proxy intercepts HTTPS traffic and presents a different and incompletely configured SSL certificate.

    To inspect the certificate chain you actually receive for `data.pypsa.org`, run:

    ```console
    pixi run openssl s_client -connect data.pypsa.org:443 -servername data.pypsa.org -showcerts </dev/null 2>/dev/null \
    | awk '/BEGIN CERTIFICATE/{i++} i{print > "/tmp/chain-cert-" i ".pem"} END{for (j=1;j<=i;j++) print "/tmp/chain-cert-" j ".pem"}' \
    | while read -r cert; do
            echo "=== $cert ==="
            openssl x509 -in "$cert" -noout -subject -issuer
            openssl x509 -in "$cert" -noout -text \
                | awk '
                    /Version:/ {print}
                    /Authority Key Identifier/ {print; getline; print}
                    /Subject Key Identifier/ {print; getline; print}
                '
        done
    ```

    This certificate should be for `data.pypsa.org` and issued by a trusted certificate authority.
    If it is not and e.g. self-signed or issued by your organisation, then this can cause the error above.

    If the certificate chain is replaced by a proxy certificate, common fixes are:

    1. Run PyPSA-Eur from a network without that proxy (for example from home).
    2. Ask your IT department for help, they can either provide a proxy bypass or ensure that their proxy certificate is correctly configured and trusted by your system.

    This is usually faster than trying to disable certificate verification or hot fixes.

??? note "I'm having trouble installing PyPSA-Eur or getting started. Where should I start?"
    The most common installation issues involve Python environment setup and solver configuration. We recommend using `pixi` for environment management. For solver setup, HiGHS is included by default for testing, but commercial solvers are supported as well. See [Installation](installation.md) for detailed platform-specific instructions, solver configuration guidance, and the legacy `conda` setup path.

??? note "My workflow is failing or producing unexpected results. How do I troubleshoot?"
    Start by running `pixi run snakemake -n` (dry-run) to validate workflow structure without execution. Then, check log files in `logs/` and verify intermediate results at each workflow stage. For persistent issues, see [Support](support.md) for community assistance channels.

??? note "I would like to develop and maintain a custom fork of PyPSA-Eur with additional features. How can I best manage this?"
    We recommend forking the PyPSA-Eur repository on GitHub and creating a dedicated branch for your custom features. Regularly sync your fork with the upstream repository to incorporate updates and improvements. Use feature branches for individual changes and consider submitting pull requests to contribute back to the main project if applicable.

??? note "What is the difference between snapshots and planning horizons?"
    **Snapshots** define the temporal resolution of system operation. They represent the individual timesteps (e.g. hourly, 4-hourly, or segments of individual lengths) over which the dispatch of technologies is optimised.

    **Planning horizons**, used for example in overnight or myopic optimisation workflows, define distinct investment periods. They either represent an individual year (e.g. 2030) or divide the long-term investment period into multiple stages (e.g. 2030, 2040, and 2050). For myopic optimisations, investments made in one horizon influence the next.

    In short: snapshots handle *operational time*, planning horizons handle *investment time*.

## Solving

??? note "My model takes much longer to solve after adding or modifying X. Why?"
    A noticeable increase in solve time usually means that the optimisation problem has become more complex. This may occur when additional variables or constraints restrict the feasible space, when the bus clustering introduces bottlenecks (for example, very small regions with high demand), or when the temporal or spatial resolution is set too high.

    For reference, PyPSA-Eur typically solves efficiently with 128 spatial clusters and 4-hour temporal resolution, assuming no custom constraints are added. If you increase either the temporal or spatial resolution, you may need to decrease the other to keep solution times manageable.

??? note "What can I do if my dispatch model reports 'infeasible or unbounded'?"
    This message typically indicates that the model cannot satisfy demand in all timesteps or that some generators have lower bounds exceeding the available load. If not already enabled, activate `load_shedding` and curtailment mode in the [solving](solving.md) configuration. These options introduce a high-cost load-shedding generator and a curtailment generator, allowing the optimisation to remain feasible.

    If enabling these options resolves the issue, examine the resulting dispatch to determine the bus regions and snapshots where the system was unable to meet demand. A simple way to quantify the total load shedding is:

    ```python
    load_shedding_i = n.generators[n.generators.carrier == "load"].index
    load_shedding = (
        n.generators_t.p[load_shedding_i]
        .mul(n.snapshot_weightings.generators, axis=0)
        .sum()
    )
    print(load_shedding)
    ```

??? note "Why does `solve_network` show the warning that the `EU` bus has no attached components?"
    The `EU` bus is introduced as an aggregate carrier bus used to track system-wide quantities—such as total fuel use (e.g. hydrogen, biomass, CO₂) across all regions. This is particularly useful when treating the respective fuel in a copperplated configuration. Because this bus is not intended to represent a physical location and no components are directly attached to it in the optimisation problem, PyPSA emits a warning. The warning is harmless and can safely be ignored.

## Model configuration & customisation

??? note "How can I change the spatial resolution of my model?"
    Several configuration options control the spatial resolution in PyPSA-Eur. Available mechanisms:

    - Adjust the number of regions via the `clusters` wildcard.
    - Restrict the geographical scope using the `countries` configuration.
      Note that the number of clusters must exceed the number of countries.
    - Configure how regions are created through the `clustering` section,
      including NUTS-based clustering (`administrative`) or custom weightings
      (`focus_weights`).
    - Provide a fully custom spatial clustering by setting `mode: custom_busshapes` and
      placing your geometry file at
      `data/busshapes/base_s_{w.clusters}_{base_network}.geojson`.

??? note "Why do some countries (e.g. GB or DK) have at least two nodes even when I specify only one cluster per country?"
    Some countries span multiple synchronous zones. Examples include Denmark (DK1 and DK2) and Great Britain versus Northern Ireland. To represent these electrical boundaries correctly, PyPSA-Eur assigns one bus region per synchronous zone, even if the configuration requests only a single cluster for the whole country.

    If you prefer a purely administrative grouping without splitting along synchronous zone boundaries, set the `clustering` option to use the NUTS classification (`administrative`) in the configuration.

??? note "How can I ensure that my custom changes are correctly implemented in the model?"
    For debugging and tracing, it is good practice to implement logs (Python's `logging` module) in your custom scripts. This allows you to track the flow of execution and identify where things might be going wrong.

    Possible steps: Use `n.statistics.energy_balance` ([docs](https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.energy_balance)) to ensure that the technology you added is expanded and dispatched. Use `n.statistics.opex` ([docs](https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.opex)) and `n.statistics.capex` ([docs](https://docs.pypsa.org/latest/api/networks/statistics/#pypsa.Network.statistics.capex)) to ensure that the overall system costs stay in a reasonable range compared to the previous (or counterfactual) run.

??? note "How can I resample an hourly time series to match my model's time resolution?"
    Use PyPSA's built-in [resampling functions](https://docs.pypsa.org/latest/api/other/common/#pypsa.common.resample_timeseries), which correctly aggregate time series and preserve all required metadata. For example, to convert an hourly profile to 4-hour resolution:

    ```python
    from pypsa.common import resample_timeseries

    ts_resampled = resample_timeseries(df, freq="4H")
    ```

    Please be aware that the dataframe's index must be a `DateTimeIndex`.

## Results & postprocessing

??? note "How can I start analysing my results (e.g. energy balances) after a successful model run?"
    PyPSA-Eur already generates a rich set of default outputs. These are stored in the `results/` directory and include maps, plots, and summary tables that provide a first overview of capacities, dispatch patterns, prices, and line loadings. Reviewing these default figures is often the quickest way to understand the high-level behaviour of your scenario.

    For more detailed analysis, use PyPSA's [built-in statistics API](https://docs.pypsa.org/latest/api/networks/statistics/), which offers a variety of functions to extract and summarise key performance indicators. These statistics can be further processed using standard Python libraries such as Pandas.

    ```python
    eb = n.statistics.energy_balance()
    print(eb)
    ```

    You can filter statistics by carrier, component, bus, or region using PyPSA's flexible filtering API. See the [documentation](https://docs.pypsa.org/stable/user-guide/statistics/#filtering) for more details.

    Other useful entry points include:

    ```python
    n.statistics.capex()      # investment costs
    n.statistics.opex()       # operational costs
    ```

    These tools help you go beyond the default outputs and systematically explore the solved system.

??? note "I observe a very high marginal electricity price in a single timestep. Is this normal?"
    This can be normal, but it is not always expected. In capacity expansion models with perfectly inelastic demand (the PyPSA-Eur default), a very high marginal price can appear when the system is operating at its limits. In such a timestep, supplying an additional unit of demand would require violating a constraint (for example, insufficient generation, storage, or transmission) or would trigger load shedding. The corresponding shadow price of the demand balance constraint then becomes large.

    However, unusually high prices can also signal modelling issues such as missing flexibility options, too few clusters causing artificial bottlenecks, an overly restrictive constraint, or incorrectly configured capacity limits. For an illustration of why marginal prices can spike when demand is perfectly inelastic, see the [demand-elasticity example](https://docs.pypsa.org/latest/examples/demand-elasticity/#perfectly-inelastic-demand-up-to-voll) in the PyPSA documentation.

## Contribution & support

??? note "How do I compile the documentation locally?"
    The documentation is built with [MkDocs](https://www.mkdocs.org/) using the dedicated `doc` Pixi environment. From the project root, run:

    ```console
    pixi run build-docs site
    ```

    This creates the HTML output in `site/html`. Open `site/html/index.html` in your browser to view the result.

??? note "How can I contribute to PyPSA-Eur?"
    We strongly welcome contributions! You can file issues or make pull requests on [GitHub](https://github.com/PyPSA/pypsa-eur). Please also refer to the [Contributing](contributing.md) section.

??? note "Where can I get help if I encounter issues?"
    Please refer to the [Support](support.md) page for various ways to reach out to us and the community, including Discord, mailing lists, and issue trackers.

??? note "Where can I report bugs or request features?"
    For bugs and feature requests, please use the repository's [GitHub issues](https://github.com/PyPSA/pypsa-eur/issues).
