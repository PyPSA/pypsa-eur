# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Snakemake PyPSA logger plugin.

Adds a logger plugin for Snakemake that provides detailed,
human-readable output tailored for PyPSA-Eur workflows.

"""

import sys
import textwrap
from dataclasses import dataclass

from snakemake.logging import (
    ColorizingTextHandler,
    DefaultFormatter,
    timestamp,
)
from snakemake_interface_logger_plugins.base import LogHandlerBase
from snakemake_interface_logger_plugins.settings import LogHandlerSettingsBase


@dataclass
class LogHandlerSettings(LogHandlerSettingsBase):
    pass


class PypsaFormatter(DefaultFormatter):
    BOLD_SEQ = "\033[1m"
    NON_BOLD_SEQ = "\033[21m"

    def format_job_info(self, msg):
        """
        Format for job_info log.

        Overrides default implementation for formatting job info messages
        to provide a more descriptive output tailored for PyPSA-Eur workflows.
        """
        output = []

        output.append(timestamp())
        if msg["rule_msg"]:
            output.append(msg["rule_msg"])
        output.append("\n".join(self._format_job_info(msg)))

        if msg.get("indent", False):
            return textwrap.indent("\n".join(output), "    ")
        return "\n".join(output)


class LogHandler(LogHandlerBase, ColorizingTextHandler):
    def __post_init__(self) -> None:
        # initialize additional attributes
        # Do not overwrite the __init__ method as this is kept in control of the base
        # class in order to simplify the update process.
        # See https://github.com/snakemake/snakemake-interface-logger-plugins/blob/main/src/snakemake_interface_logger_plugins/base.py # noqa: E501
        # for attributes of the base class.
        # In particular, the settings of above LogHandlerSettings class are accessible via
        # self.settings.
        # You also have access to self.common_settings here, which are logging settings supplied by the caller in the form of OutputSettingsLoggerInterface. # noqa: E501
        # See https://github.com/snakemake/snakemake-interface-logger-plugins/blob/main/src/snakemake_interface_logger_plugins/settings.py for more details # noqa: E501
        ColorizingTextHandler.__init__(
            self,
            nocolor=self.common_settings.nocolor,
            stream=sys.stdout if self.common_settings.stdout else sys.stderr,
        )
        self.setFormatter(
            PypsaFormatter(
                self.common_settings.quiet, self.common_settings.show_failed_logs
            )
        )

    # Here you can override logging.Handler methods to customize logging behavior.
    # Only an implementation of the emit() method is required. See the Python logging
    # documentation for details:
    # https://docs.python.org/3/library/logging.html#handler-objects

    # LogRecords from Snakemake carry contextual information in the record's attributes
    # Of particular interest is the 'event' attribute, which indicates the type of log information contained
    # See https://github.com/snakemake/snakemake-interface-logger-plugins/blob/2ab84cb31f0b92cf0b7ee3026e15d1209729d197/src/snakemake_interface_logger_plugins/common.py#L33 # noqa: E501
    # For examples on parsing LogRecords, see https://github.com/cademirch/snakemake-logger-plugin-snkmt/blob/main/src/snakemake_logger_plugin_snkmt/parsers.py # noqa: E501

    def emit(self, record):
        # Actually emit the record. Typically this will call self.format(record) to
        # convert the record to a formatted string. The result could then be written to
        # a stream or file.
        ColorizingTextHandler.emit(self, record)

    @property
    def writes_to_stream(self) -> bool:
        # Whether this plugin writes to stderr/stdout.
        # If your plugin writes to stderr/stdout, return
        # true so that Snakemake disables its stderr logging.
        return True

    @property
    def writes_to_file(self) -> bool:
        # Whether this plugin writes to a file.
        # If your plugin writes log output to a file, return
        # true so that Snakemake can report your logfile path at workflow end.
        # NOTE: Handlers that return True must provide a baseFilename attribute
        # containing the file path.
        return False

    @property
    def has_filter(self) -> bool:
        # Whether this plugin attaches its own filter.
        # Return true if your plugin provides custom log filtering logic.
        # If false is returned, Snakemake's DefaultFilter will be attached see: https://github.com/snakemake/snakemake/blob/960f6a89eaa31da6014e810dfcf08f635ac03a6e/src/snakemake/logging.py#L372 # noqa: E501
        # See https://docs.python.org/3/library/logging.html#filter-objects for info on how to define and attach a Filter
        return False

    @property
    def has_formatter(self) -> bool:
        # Whether this plugin attaches its own formatter.
        # Return true if your plugin provides custom log formatting logic.
        # If false is returned, Snakemake's Defaultformatter will be attached see: https://github.com/snakemake/snakemake/blob/960f6a89eaa31da6014e810dfcf08f635ac03a6e/src/snakemake/logging.py#L132 # noqa: E501
        # See https://docs.python.org/3/library/logging.html#formatter-objects for info on how to define and attach a Formatter
        return True

    @property
    def needs_rulegraph(self) -> bool:
        # Whether this plugin requires the DAG rulegraph.
        # Return true if your plugin needs access to the workflow's
        # directed acyclic graph for logging purposes.
        return False
