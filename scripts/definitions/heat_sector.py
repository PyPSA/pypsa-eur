# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum


class HeatSector(Enum):
    """
    Enumeration class representing different heat sectors.

    Attributes:
        RESIDENTIAL (str): Represents the residential heat sector.
        SERVICES (str): Represents the services heat sector.
    """

    RESIDENTIAL = "residential"
    SERVICES = "services"

    def __str__(self) -> str:
        """
        Returns the string representation of the heat sector.

        Returns:
            str: The string representation of the heat sector.
        """
        return self.value
