# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from pathlib import Path
import yaml
from os.path import normpath, exists, join
from shutil import copyfile, move, rmtree
from snakemake.utils import min_version

min_version("8.11")

from scripts._helpers import (
    get_rdir,
    get_scenarios,
    get_shadow,
    path_provider,
)
from scripts.lib.validation.config import normalize_config, validate_config


configfile: "config/config.default.yaml"
configfile: "config/plotting.default.yaml"


if Path("config/config.yaml").exists():

    configfile: "config/config.yaml"


validated = validate_config(config)
normalize_config(config, validated)

run = config["run"]
scenarios = get_scenarios(run)
RDIR = get_rdir(run)
shadow_config = get_shadow(run)

shared_resources = run["shared_resources"]["policy"]
exclude_from_shared = run["shared_resources"]["exclude"]
logs = path_provider("logs/", RDIR, shared_resources, exclude_from_shared)
_benchmark_provider = path_provider(
    "benchmarks/", RDIR, shared_resources, exclude_from_shared
)
resources = path_provider("resources/", RDIR, shared_resources, exclude_from_shared)


def benchmarks(fn):
    """Return a benchmark file path, even if legacy directories already exist."""

    path = Path(_benchmark_provider(fn))
    if path.is_dir():
        # Preserve existing benchmark directories by placing the file inside.
        path = path / "benchmark.tsv"
    elif not path.suffix:
        path = path.with_suffix(".tsv")
    return str(path)


RESULTS = "results/" + RDIR


localrules:
    purge,


wildcard_constraints:
    horizon=r"[0-9]{4}",


include: "rules/common.smk"


# Data constants
OSM_DATASET = dataset_version("osm")


include: "rules/collect.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/compose.smk"
include: "rules/solve.smk"
include: "rules/postprocess.smk"
include: "rules/development.smk"


# Define output categories based on foresight mode
# This follows the same pattern as postprocess.smk for consistency

# Core outputs that always run
CORE_OUTPUTS = [
    RESULTS + "graphs/costs.svg",
    RESULTS + "graphs/energy.svg",
    RESULTS + "graphs/balances_energy.svg",
]

# Network and timeseries plots (excluded for perfect foresight)
if config["foresight"] != "perfect":
    NETWORK_PLOT_OUTPUTS = [
        resources("maps/base_network.pdf"),
        resources("maps/clustered_network.pdf"),
        RESULTS + "maps/static/power_network_{horizon}.pdf",
        RESULTS + "graphs/cop_profiles_{horizon}.html",
    ]
    TIMESERIES_OUTPUTS = [
        RESULTS + "graphs/balance_timeseries_{horizon}",
        RESULTS + "graphs/heatmap_timeseries_{horizon}",
    ]
else:
    NETWORK_PLOT_OUTPUTS = []
    TIMESERIES_OUTPUTS = []

# Myopic-specific outputs
if config["foresight"] == "myopic":
    MYOPIC_OUTPUTS = [RESULTS + "csvs/cumulative_costs.csv"]
else:
    MYOPIC_OUTPUTS = []


def get_sector_network_plots(w):
    """Returns sector-specific network plots if enabled and not perfect foresight."""
    if config["foresight"] == "perfect":
        return []

    plots = []
    if config_provider("sector", "H2_network")(w):
        plots.extend(
            expand(
                RESULTS + "maps/static/h2_network_{horizon}.pdf",
                horizon=config["planning_horizons"],
                run=config["run"]["name"],
            )
        )
    if config_provider("sector", "gas_network")(w):
        plots.extend(
            expand(
                RESULTS + "maps/static/ch4_network_{horizon}.pdf",
                horizon=config["planning_horizons"],
                run=config["run"]["name"],
            )
        )
    return plots


def get_balance_map_plots(w):
    """Returns balance map plots (static + interactive) if configured and not perfect foresight."""
    if config["foresight"] == "perfect":
        return []

    plots = []

    # Static balance maps (PDF)
    static_carriers = config_provider("plotting", "balance_map", "bus_carriers")(w)
    if static_carriers:
        plots.extend(
            expand(
                RESULTS + "maps/static/balance_map_{carrier}_{horizon}.pdf",
                horizon=config["planning_horizons"],
                run=config["run"]["name"],
                carrier=static_carriers,
            )
        )

    # Interactive balance maps (HTML)
    interactive_carriers = config_provider(
        "plotting", "balance_map_interactive", "bus_carriers"
    )(w)
    if interactive_carriers:
        plots.extend(
            expand(
                RESULTS + "maps/interactive/balance_map_{carrier}_{horizon}.html",
                horizon=config["planning_horizons"],
                run=config["run"]["name"],
                carrier=interactive_carriers,
            )
        )

    return plots


rule all:
    input:
        expand(CORE_OUTPUTS, run=config["run"]["name"]),
        (
            expand(
                NETWORK_PLOT_OUTPUTS,
                run=config["run"]["name"],
                horizon=config["planning_horizons"],
            )
            if NETWORK_PLOT_OUTPUTS
            else []
        ),
        (
            expand(
                TIMESERIES_OUTPUTS,
                run=config["run"]["name"],
                horizon=config["planning_horizons"],
            )
            if TIMESERIES_OUTPUTS
            else []
        ),
        (expand(MYOPIC_OUTPUTS, run=config["run"]["name"]) if MYOPIC_OUTPUTS else []),
        get_sector_network_plots,
        get_balance_map_plots,
    default_target: True


rule create_scenarios:
    output:
        config["run"]["scenarios"]["file"],
    script:
        "config/create_scenarios.py"


rule purge:
    run:
        import builtins

        do_purge = builtins.input(
            "Do you really want to delete all generated files?\n"
            "\t* resources\n"
            "\t* results\n"
            "\t* docs\n"
            "Downloaded files are kept.\n"
            "Delete all files in the folders above? [y/N] "
        )
        if do_purge == "y":

            # Remove the directories and recreate them with .gitkeep
            for dir_path in ["resources/", "results/"]:
                rmtree(dir_path, ignore_errors=True)
                Path(dir_path).mkdir(parents=True, exist_ok=True)
                (Path(dir_path) / ".gitkeep").touch()

            rmtree("doc/_build", ignore_errors=True)
            print(
                "Purging all generated resources, results and docs. Downloads are kept."
            )
        else:
            raise Exception(f"Input {do_purge}. Aborting purge.")


rule dump_graph_config:
    """Dump the current Snakemake configuration to a YAML file for graph generation."""
    output:
        config_file=temp(resources("dag_final_config.yaml")),
    run:
        import yaml

        with open(output.config_file, "w") as f:
            yaml.dump(config, f)


rule rulegraph:
    """Generates Rule DAG in DOT, PDF, PNG, and SVG formats using the final configuration."""
    message:
        "Creating RULEGRAPH dag in multiple formats using the final configuration."
    input:
        config_file=rules.dump_graph_config.output.config_file,
    output:
        dot=resources("dag_rulegraph.dot"),
        pdf=resources("dag_rulegraph.pdf"),
        png=resources("dag_rulegraph.png"),
        svg=resources("dag_rulegraph.svg"),
    shell:
        r"""
        # Generate DOT file using nested snakemake with the dumped final config
        echo "[Rule rulegraph] Using final config file: {input.config_file}"
        snakemake --rulegraph --configfile {input.config_file} --quiet | sed -n "/digraph/,\$p" > {output.dot}

        # Generate visualizations from the DOT file
        if [ -s {output.dot} ]; then
            dot -c

            echo "[Rule rulegraph] Generating PDF from DOT"
            dot -Tpdf -o {output.pdf} {output.dot} || {{ echo "Error: Failed to generate PDF. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule rulegraph] Generating PNG from DOT"
            dot -Tpng -o {output.png} {output.dot} || {{ echo "Error: Failed to generate PNG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule rulegraph] Generating SVG from DOT"
            dot -Tsvg -o {output.svg} {output.dot} || {{ echo "Error: Failed to generate SVG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule rulegraph] Successfully generated all formats."
        else
            echo "[Rule rulegraph] Error: Failed to generate valid DOT content." >&2
            exit 1
        fi
        """


rule filegraph:
    """Generates File DAG in DOT, PDF, PNG, and SVG formats using the final configuration."""
    message:
        "Creating FILEGRAPH dag in multiple formats using the final configuration."
    input:
        config_file=rules.dump_graph_config.output.config_file,
    output:
        dot=resources("dag_filegraph.dot"),
        pdf=resources("dag_filegraph.pdf"),
        png=resources("dag_filegraph.png"),
        svg=resources("dag_filegraph.svg"),
    shell:
        r"""
        # Generate DOT file using nested snakemake with the dumped final config
        echo "[Rule filegraph] Using final config file: {input.config_file}"
        snakemake --filegraph all --configfile {input.config_file} --quiet | sed -n "/digraph/,\$p" > {output.dot}

        # Generate visualizations from the DOT file
        if [ -s {output.dot} ]; then
            echo "[Rule filegraph] Generating PDF from DOT"
            dot -Tpdf -o {output.pdf} {output.dot} || {{ echo "Error: Failed to generate PDF. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule filegraph] Generating PNG from DOT"
            dot -Tpng -o {output.png} {output.dot} || {{ echo "Error: Failed to generate PNG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule filegraph] Generating SVG from DOT"
            dot -Tsvg -o {output.svg} {output.dot} || {{ echo "Error: Failed to generate SVG. Is graphviz installed?" >&2; exit 1; }}

            echo "[Rule filegraph] Successfully generated all formats."
        else
            echo "[Rule filegraph] Error: Failed to generate valid DOT content." >&2
            exit 1
        fi
        """


rule doc:
    message:
        "Build documentation."
    output:
        directory("doc/_build"),
    shell:
        "pixi run build-docs {output} html"


rule sync:
    params:
        cluster=f"{config['remote']['ssh']}:{config['remote']['path']}",
    shell:
        """
        rsync -uvarh --ignore-missing-args --files-from=.sync-send . {params.cluster}
        rsync -uvarh --no-g {params.cluster}/resources . || echo "No resources directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/results . || echo "No results directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/logs . || echo "No logs directory, skipping rsync"
        """


rule sync_dry:
    params:
        cluster=f"{config['remote']['ssh']}:{config['remote']['path']}",
    shell:
        """
        rsync -uvarh --ignore-missing-args --files-from=.sync-send . {params.cluster} -n
        rsync -uvarh --no-g {params.cluster}/resources . -n || echo "No resources directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/results . -n || echo "No results directory, skipping rsync"
        rsync -uvarh --no-g {params.cluster}/logs . -n || echo "No logs directory, skipping rsync"
        """
