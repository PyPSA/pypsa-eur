
import logging
logger = logging.getLogger(__name__)

import pandas as pd
from pathlib import Path

def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import snakemake as sm
    import os
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert Path.cwd().resolve() == script_dir, \
      f'mock_snakemake has to be run from the repository scripts directory {script_dir}'
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    workflow = sm.Workflow(snakefile)
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(job.input, job.output, job.params, job.wildcards,
                          job.threads, job.resources, job.log,
                          job.dag.workflow.config, job.rule.name, None,)
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)
    return snakemake

#https://stackoverflow.com/questions/20833344/fix-invalid-polygon-in-shapely
#https://stackoverflow.com/questions/13062334/polygon-intersection-error-in-shapely-shapely-geos-topologicalerror-the-opera
#https://shapely.readthedocs.io/en/latest/manual.html#object.buffer
def clean_invalid_geometries(geometries):
    """Fix self-touching or self-crossing polygons; these seem to appear
due to numerical problems from writing and reading, since the geometries
are valid before being written in pypsa-eur/scripts/cluster_network.py"""
    for i,p in geometries.items():
        if not p.is_valid:
            logger.warning(f'Clustered region {i} had an invalid geometry, fixing using zero buffer.')
            geometries[i] = p.buffer(0)
