# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Overpass API configuration.
"""

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _UserAgentConfig(ConfigModel):
    """Configuration for `overpass_api.user_agent` settings."""

    project_name: str = Field(
        "PyPSA-Eur",
        description="Project name used to identify the user agent of the Overpass API requests.",
    )
    email: str = Field(
        "contact@pypsa.org",
        description="Contact email address for the project using the Overpass API.",
    )
    website: str = Field(
        "https://github.com/PyPSA/pypsa-eur",
        description="Website URL for the project using the Overpass API.",
    )


class OverpassApiConfig(BaseModel):
    """Configuration for `overpass_api` settings."""

    url: str = Field(
        "https://overpass-api.de/api/interpreter",
        description="Overpass API endpoint URL. See `Overpass API Wiki <https://wiki.openstreetmap.org/wiki/Overpass_API#Public_Overpass_API_instances>`_ for available public instances.",
    )
    max_tries: int = Field(
        5,
        description="Maximum retry attempts for Overpass API requests. Please be respectful to the Overpass API fair use policy of the individual instances.",
    )
    timeout: int = Field(
        600,
        description="Timeout in seconds for Overpass API requests.",
    )
    user_agent: _UserAgentConfig = Field(
        default_factory=_UserAgentConfig,
        description="Please provide your own user agent details when using the Overpass API,so the instance operators can contact you if needed.",
    )
