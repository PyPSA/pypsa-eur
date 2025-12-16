# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Run configuration block.

See # docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#run
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class _ScenariosConfig(ConfigModel):
    """Configuration for `run.scenarios` level."""

    enable: bool = Field(
        False,
        description="Switch to select whether workflow should generate scenarios based on ``file``.",
    )
    file: str = Field(
        "config/scenarios.yaml",
        description="Path to the scenario yaml file. The scenario file contains config overrides for each scenario. In order to be taken account, ``run: scenarios`` has to be set to ``true`` and ``run: name`` has to be a subset of top level keys given in the scenario file. In order to automatically create a `scenario.yaml` file based on a combination of settings, alter and use the ``config/create_scenarios.py`` script in the ``config`` directory.",
        examples=["config/scenarios.yaml"],
    )


class _SharedResourcesConfig(ConfigModel):
    """Configuration for `run.shared_resources` level."""

    policy: bool | str = Field(
        False,
        description="Boolean switch to select whether resources should be shared across runs. If a string is passed, this is used as a subdirectory name for shared resources. If set to 'base', only resources before creating the elec.nc file are shared.",
        examples=[False, "base"],
    )
    exclude: list[str] = Field(
        default_factory=list,
        description="For the case shared_resources=base, specify additional files that should not be shared across runs.",
    )


class RunConfig(ConfigModel):
    """Configuration for top level `run` settings."""

    prefix: str = Field(
        "",
        description="Prefix for the run name which is used as a top-layer directory name in the results and resources folders.",
    )
    name: str | list[str] = Field(
        "",
        description="Specify a name for your run. Results will be stored under this name. If ``scenario: enable:`` is set to ``true``, the name must contain a subset of scenario names defined in ``scenario: file:``. If the name is 'all', all defined scenarios will be run.",
        examples=["my-pypsa-eur-run"],
    )

    scenarios: _ScenariosConfig = Field(
        default_factory=_ScenariosConfig,
        description="Configuration for running multiple scenarios",
    )

    disable_progressbar: bool = Field(
        False, description="Switch to select whether progressbar should be disabled."
    )

    shared_resources: _SharedResourcesConfig = Field(
        default_factory=_SharedResourcesConfig,
        description="Shared resources configuration for parallel execution",
    )

    use_shadow_directory: bool = Field(
        False,
        description="Set to ``true`` (default) if snakemake shadow directories (``shallow``) should be used. Set to ``false`` if problems occur.",
        examples=[True],
    )
