# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum

### NOTIZEN: Analog hierzu eine definition für heating_components ermöglichen und dabei einbauen, dass if config true nehmen wir für die Komponeten an, dass sie als booster fungieren können
# dann loopen wir in prepare sector newtork über die komponenten

class HeatSystemType(Enum):
    """
    Enumeration representing different types of heat systems.
    """

    URBAN_CENTRAL = "urban central"
    URBAN_DECENTRAL = "urban decentral"
    RURAL = "rural"

    def __str__(self) -> str:
        """
        Returns the string representation of the heat system type.

        Returns:
            str: The string representation of the heat system type.
        """
        return self.value

    @property
    def is_central(self) -> bool:
        """
        Returns whether the heat system type is central.

        Returns:
            bool: True if the heat system type is central, False otherwise.
        """
        return self == HeatSystemType.URBAN_CENTRAL
