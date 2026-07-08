# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Stochastic scenarios configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#stochastic-scenarios
"""

from pathlib import Path

from pydantic import Field, field_validator, model_validator

from scripts.lib.validation.config._base import ConfigModel


class StochasticScenariosExportConfig(ConfigModel):
    """Configuration for exporting deterministic views of stochastic solutions."""

    average: bool = Field(
        True,
        description=(
            "Export a probability-weighted average deterministic network view from "
            "the stochastic solution."
        ),
    )

    scenarios: bool = Field(
        False,
        description=(
            "Export one deterministic network view for each stochastic scenario."
        ),
    )


class StochasticScenariosPostprocessConfig(ConfigModel):
    """Configuration for post-processing stochastic solutions."""

    use_average: bool = Field(
        True,
        description=(
            "Use the probability-weighted average deterministic network view as "
            "input for standard post-processing rules when stochastic scenarios "
            "are enabled."
        ),
    )

    scenarios: bool = Field(
        False,
        description=(
            "Run scenario-specific post-processing for all stochastic scenarios."
        ),
    )


class StochasticScenariosConfig(ConfigModel):
    """Configuration for stochastic scenarios."""

    enable: bool = Field(
        False,
        description="Whether to build and solve a stochastic PyPSA network.",
    )

    file: str = Field(
        "config/stochastic_scenarios.yaml",
        description="Path to the external YAML file defining stochastic scenarios.",
    )

    export: StochasticScenariosExportConfig = Field(
        default_factory=StochasticScenariosExportConfig,
        description="Export options for stochastic solutions.",
    )

    postprocess: StochasticScenariosPostprocessConfig = Field(
        default_factory=StochasticScenariosPostprocessConfig,
        description="Post-processing options for stochastic solutions.",
    )

    @field_validator("file")
    @classmethod
    def file_must_be_yaml(cls, value: str) -> str:
        """Check that the stochastic scenario file has a YAML extension."""
        if not value.endswith(".yaml") and not value.endswith(".yml"):
            raise ValueError("stochastic_scenarios.file must be a YAML file.")
        return value

    @model_validator(mode="after")
    def check_average_export_for_postprocess(self):
        """Ensure average post-processing has the required exported network."""
        if self.enable and self.postprocess.use_average and not self.export.average:
            raise ValueError(
                "stochastic_scenarios.postprocess.use_average requires "
                "stochastic_scenarios.export.average to be true."
            )
        return self

    @model_validator(mode="after")
    def check_scenario_export_for_postprocess(self):
        """Ensure scenario-specific post-processing has scenario exports enabled."""
        if self.enable and self.postprocess.scenarios and not self.export.scenarios:
            raise ValueError(
                "stochastic_scenarios.postprocess.scenarios requires "
                "stochastic_scenarios.export.scenarios to be true."
            )
        return self