# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Foresight configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#foresight
"""

from typing import Literal

from pydantic import Field, RootModel


class ForesightConfig(RootModel[Literal["overnight", "myopic", "perfect"]]):
    """Configuration for `foresight` settings."""

    root: Literal["overnight", "myopic", "perfect"] = Field(
        "overnight",
        description="See Foresight Options for detail explanations.",
    )
