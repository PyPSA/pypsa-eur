# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Transformers configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#transformers
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class TransformersConfig(ConfigModel):
    """Configuration for `transformers` settings."""

    x: float = Field(
        0.1,
        description="Series reactance in per unit (p.u.), using `s_nom` as base power of the transformer. Overwritten if `type` is specified.",
    )
    s_nom: float = Field(
        2000.0,
        description="Limit of apparent power which can pass through branch (MVA). Overwritten if `type` is specified.",
    )
    type: str = Field(
        "",
        description="Specifies transformer types to assume for the transformers of the ENTSO-E grid extraction.",
    )
