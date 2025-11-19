# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Top-level configuration fields.

See https://pypsa-eur.readthedocs.io/en/latest/configuration.html#top-level-configuration.

"""

from typing import Literal

from pydantic import BaseModel, Field


class LoggingConfig(BaseModel):
    """Configuration for top level `logging` settings."""

    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = Field(
        "INFO",
        description="Restrict console outputs to all infos, warning or errors only",
    )
    format: str = Field(
        "%(levelname)s:%(name)s:%(message)s",
        description="Custom format for log messages. See `LogRecord <https://docs.python.org/3/library/logging.html#logging.LogRecord>`_ attributes.",
    )


class RemoteConfig(BaseModel):
    """Configuration for top level `remote` settings."""

    ssh: str = Field(
        "",
        description="Optionally specify the SSH of a remote cluster to be synchronized.",
    )
    path: str = Field(
        "",
        description="Optionally specify the file path within the remote cluster to be synchronized.",
    )


class TopLevelConfig(BaseModel):
    """Top level configuration."""

    version: str = Field(
        "v2025.07.0", description="Version of PyPSA-Eur. Descriptive only."
    )

    tutorial: bool = Field(
        False,
        description="Switch to retrieve the tutorial data set instead of the full data set.",
    )

    logging: LoggingConfig = Field(
        default_factory=LoggingConfig,
        description="Logging configuration for the workflow",
    )

    remote: RemoteConfig = Field(
        default_factory=RemoteConfig,
        description="Configuration for remote workflow execution",
    )
