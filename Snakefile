# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from pathlib import Path
import yaml
from os.path import normpath, exists, join
from shutil import copyfile, move, rmtree
from dotenv import load_dotenv
from snakemake.utils import min_version

load_dotenv()

min_version("8.11")

from scripts._helpers import (
    get_rdir,
    get_scenarios,
    get_shadow,
    path_provider,
    script_path_provider,
)
from scripts.lib.validation.config import validate_config


configfile: "config/config.default.yaml"
configfile: "config/plotting.default.yaml"


if Path("config/config.yaml").exists():

    configfile: "config/config.yaml"


validate_config(config)

run = config["run"]
scenarios = get_scenarios(run)
RDIR = get_rdir(run)
shadow_config = get_shadow(run)

shared_resources = run["shared_resources"]["policy"]
exclude_from_shared = run["shared_resources"]["exclude"]
logs = path_provider("logs/", RDIR, shared_resources, exclude_from_shared)
benchmarks = path_provider("benchmarks/", RDIR, shared_resources, exclude_from_shared)
resources = path_provider("resources/", RDIR, shared_resources, exclude_from_shared)
scripts = script_path_provider(Path(workflow.snakefile).parent)

RESULTS = "results/" + RDIR


localrules:
    purge,


wildcard_constraints:
    clusters="[0-9]+(m|c)?|all|adm",
    opts=r"[-+a-zA-Z0-9\.]*",
    sector_opts=r"[-+a-zA-Z0-9\.\s]*",
    planning_horizons=r"[0-9]{4}",


include: "rules/common.smk"


# Data constants
OSM_DATASET = dataset_version("osm")


include: "rules/collect.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/solve_electricity.smk"
include: "rules/postprocess.smk"
include: "rules/development.smk"


if config["foresight"] == "overnight":

    include: "rules/solve_overnight.smk"


if config["foresight"] == "myopic":

    include: "rules/solve_myopic.smk"


if config["foresight"] == "perfect":

    include: "rules/solve_perfect.smk"


rule all:
    input:
        expand(RESULTS + "graphs/costs.pdf", run=config["run"]["name"]),
        expand(resources("maps/power-network.pdf"), run=config["run"]["name"]),
        expand(
            resources("maps/power-network-s-{clusters}.pdf"),
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        # COP profiles plots
        expand(
            RESULTS + "graphs/cop_profiles_s_{clusters}_{planning_horizons}.html",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-h2_network_{planning_horizons}.pdf"
                if config_provider("sector", "H2_network")(w)
                else []
            ),
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}-ch4_network_{planning_horizons}.pdf"
                if config_provider("sector", "gas_network")(w)
                else []
            ),
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: expand(
            (
                RESULTS + "csvs/cumulative_costs.csv"
                if config_provider("foresight")(w) == "myopic"
                else []
            ),
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "graphics/balance_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "graphics/heatmap_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        # Explicitly list heat source types for temperature maps
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_temperature_map_river_water.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "river_water"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_temperature_map_sea_water.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "sea_water"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_temperature_map_ambient_air.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "air"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        # Only river_water has energy maps
        lambda w: expand(
            (
                RESULTS
                + "maps/static/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}-heat_source_energy_map_river_water.html"
                if config_provider("plotting", "enable_heat_source_maps")(w)
                and "river_water"
                in config_provider("sector", "heat_pump_sources", "urban central")(w)
                else []
            ),
            **config["scenario"],
            run=config["run"]["name"],
        ),
        expand(
            RESULTS
            + "graphics/balance_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "graphics/heatmap_timeseries/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        expand(
            RESULTS
            + "graphics/interactive_bus_balance/s_{clusters}_{opts}_{sector_opts}_{planning_horizons}",
            run=config["run"]["name"],
            **config["scenario"],
        ),
        lambda w: balance_map_paths("static", w),
        lambda w: balance_map_paths("interactive", w),
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
