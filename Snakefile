# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from pathlib import Path
import yaml
from os.path import normpath, exists
from shutil import copyfile, move, rmtree
from snakemake.utils import min_version

min_version("8.11")

from scripts._helpers import path_provider, copy_default_files, get_scenarios, get_rdir


copy_default_files(workflow)


configfile: "config/config.default.yaml"
configfile: "config/config.yaml"


run = config["run"]
scenarios = get_scenarios(run)
RDIR = get_rdir(run)

shared_resources = run["shared_resources"]["policy"]
exclude_from_shared = run["shared_resources"]["exclude"]
logs = path_provider("logs/", RDIR, shared_resources, exclude_from_shared)
benchmarks = path_provider("benchmarks/", RDIR, shared_resources, exclude_from_shared)
resources = path_provider("resources/", RDIR, shared_resources, exclude_from_shared)

CDIR = "" if run["shared_cutouts"] else RDIR
RESULTS = "results/" + RDIR


localrules:
    purge,


wildcard_constraints:
    clusters="[0-9]+(m|c)?|all",
    ll=r"(v|c)([0-9\.]+|opt)",
    opts=r"[-+a-zA-Z0-9\.]*",
    sector_opts=r"[-+a-zA-Z0-9\.\s]*",


include: "rules/common.smk"
include: "rules/collect.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/solve_electricity.smk"
include: "rules/postprocess.smk"
include: "rules/validate.smk"
include: "rules/development.smk"


if config["foresight"] == "overnight":

    include: "rules/solve_overnight.smk"


if config["foresight"] == "myopic":

    include: "rules/solve_myopic.smk"


if config["foresight"] == "perfect":

    include: "rules/solve_perfect.smk"


rule all:
    input:
        expand(RESULTS + "graphs/costs.svg", run=config["run"]["name"]),
    default_target: True


rule create_scenarios:
    output:
        config["run"]["scenarios"]["file"],
    conda:
        "envs/retrieve.yaml"
    script:
        "config/create_scenarios.py"


rule purge:
    run:
        import builtins

        do_purge = builtins.input(
            "Do you really want to delete all generated resources, \nresults and docs (downloads are kept)? [y/N] "
        )
        if do_purge == "y":
            rmtree("resources/", ignore_errors=True)
            rmtree("results/", ignore_errors=True)
            rmtree("doc/_build", ignore_errors=True)
            print("Purging generated resources, results and docs. Downloads are kept.")
        else:
            raise Exception(f"Input {do_purge}. Aborting purge.")


rule dag:
    message:
        "Creating DAG of workflow."
    output:
        dot=resources("dag.dot"),
        pdf=resources("dag.pdf"),
        png=resources("dag.png"),
    conda:
        "envs/environment.yaml"
    shell:
        r"""
        snakemake --rulegraph all | sed -n "/digraph/,\$p" > {output.dot}
        dot -Tpdf -o {output.pdf} {output.dot}
        dot -Tpng -o {output.png} {output.dot}
        """


rule doc:
    message:
        "Build documentation."
    output:
        directory("doc/_build"),
    shell:
        "make -C doc html"


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
