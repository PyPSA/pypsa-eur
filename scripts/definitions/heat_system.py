# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum

from definitions.heat_sector import HeatSector
from definitions.heat_system_type import HeatSystemType


class HeatSystem(Enum):
    """
    Enumeration representing different heat systems.

    Attributes
    ----------
    RESIDENTIAL_RURAL : str
        Heat system for residential areas in rural locations.
    SERVICES_RURAL : str
        Heat system for service areas in rural locations.
    RESIDENTIAL_URBAN_DECENTRAL : str
        Heat system for residential areas in urban decentralized locations.
    SERVICES_URBAN_DECENTRAL : str
        Heat system for service areas in urban decentralized locations.
    URBAN_CENTRAL : str
        Heat system for urban central areas.

    Methods
    -------
    __str__()
        Returns the string representation of the heat system.
    central_or_decentral()
        Returns whether the heat system is central or decentralized.
    system_type()
        Returns the type of the heat system.
    sector()
        Returns the sector of the heat system.
    rural()
        Returns whether the heat system is for rural areas.
    urban_decentral()
        Returns whether the heat system is for urban decentralized areas.
    urban()
        Returns whether the heat system is for urban areas.
    heat_demand_weighting(urban_fraction=None, dist_fraction=None)
        Calculates the heat demand weighting based on urban fraction and distribution fraction.
    heat_pump_costs_name(heat_source)
        Generates the name for the heat pump costs based on the heat source.
    """

    RESIDENTIAL_RURAL = "residential rural"
    SERVICES_RURAL = "services rural"
    RESIDENTIAL_URBAN_DECENTRAL = "residential urban decentral"
    SERVICES_URBAN_DECENTRAL = "services urban decentral"
    URBAN_CENTRAL = "urban central"

    def __init__(self, *args):
        super().__init__(*args)

    def __str__(self) -> str:
        """
        Returns the string representation of the heat system.

        Returns
        -------
        str
            The string representation of the heat system.
        """
        return self.value

    @property
    def central_or_decentral(self) -> str:
        """
        Returns whether the heat system is central or decentralized.

        Returns
        -------
        str
            "central" if the heat system is central, "decentral" otherwise.
        """
        if self == HeatSystem.URBAN_CENTRAL:
            return "central"
        else:
            return "decentral"

    @property
    def system_type(self) -> HeatSystemType:
        """
        Returns the type of the heat system.

        Returns
        -------
        str
            The type of the heat system.

        Raises
        ------
        RuntimeError
            If the heat system is invalid.
        """
        if self == HeatSystem.URBAN_CENTRAL:
            return HeatSystemType.URBAN_CENTRAL
        elif (
            self == HeatSystem.RESIDENTIAL_URBAN_DECENTRAL
            or self == HeatSystem.SERVICES_URBAN_DECENTRAL
        ):
            return HeatSystemType.URBAN_DECENTRAL
        elif self == HeatSystem.RESIDENTIAL_RURAL or self == HeatSystem.SERVICES_RURAL:
            return HeatSystemType.RURAL
        else:
            raise RuntimeError(f"Invalid heat system: {self}")

    @property
    def sector(self) -> HeatSector:
        """
        Returns the sector of the heat system.

        Returns
        -------
        HeatSector
            The sector of the heat system.
        """
        if (
            self == HeatSystem.RESIDENTIAL_RURAL
            or self == HeatSystem.RESIDENTIAL_URBAN_DECENTRAL
        ):
            return HeatSector.RESIDENTIAL
        elif (
            self == HeatSystem.SERVICES_RURAL
            or self == HeatSystem.SERVICES_URBAN_DECENTRAL
        ):
            return HeatSector.SERVICES
        else:
            "tot"

    @property
    def is_rural(self) -> bool:
        """
        Returns whether the heat system is for rural areas.

        Returns
        -------
        bool
            True if the heat system is for rural areas, False otherwise.
        """
        if self == HeatSystem.RESIDENTIAL_RURAL or self == HeatSystem.SERVICES_RURAL:
            return True
        else:
            return False

    @property
    def is_urban_decentral(self) -> bool:
        """
        Returns whether the heat system is for urban decentralized areas.

        Returns
        -------
        bool
            True if the heat system is for urban decentralized areas, False otherwise.
        """
        if (
            self == HeatSystem.RESIDENTIAL_URBAN_DECENTRAL
            or self == HeatSystem.SERVICES_URBAN_DECENTRAL
        ):
            return True
        else:
            return False

    @property
    def is_urban(self) -> bool:
        """
        Returns whether the heat system is for urban areas.

        Returns
        -------
        bool True if the heat system is for urban areas, False otherwise.
        """
        return not self.is_rural

    def heat_demand_weighting(self, urban_fraction=None, dist_fraction=None) -> float:
        """
        Calculates the heat demand weighting based on urban fraction and
        distribution fraction.

        Parameters
        ----------
        urban_fraction : float, optional
            The fraction of urban heat demand.
        dist_fraction : float, optional
            The fraction of distributed heat demand.

        Returns
        -------
        float
            The heat demand weighting.

        Raises
        ------
        RuntimeError
            If the heat system is invalid.
        """
        if "rural" in self.value:
            return 1 - urban_fraction
        elif "urban central" in self.value:
            return dist_fraction
        elif "urban decentral" in self.value:
            return urban_fraction - dist_fraction
        else:
            raise RuntimeError(f"Invalid heat system: {self}")

    def heat_pump_costs_name(self, heat_source: str) -> str:
        """
        Generates the name for the heat pump costs based on the heat source and
        system.
        Used to retrieve data from `technology-data <https://github.com/PyPSA/technology-data>`.

        Parameters
        ----------
        heat_source : str
            The heat source.

        Returns
        -------
        str
            The name for the heat pump costs.
        """
        return f"{self.central_or_decentral} {heat_source}-sourced heat pump"

    def heat_source_costs_name(self, heat_source: str) -> str:
        """
        Generates the name for direct source utilisation costs based on the heat source and
        system.
        Used to retrieve data from `technology-data <https://github.com/PyPSA/technology-data>`.

        Parameters
        ----------
        heat_source : str
            The heat source.

        Returns
        -------
        str
            The name for the technology-data costs.
        """
        return f"{self.central_or_decentral} {heat_source} heat source"

    @property
    def resistive_heater_costs_name(self) -> str:
        """
        Generates the name for the resistive heater costs based on the heat
        system.
        Used to retrieve data from `technology-data <https://github.com/PyPSA/technology-data>`.

        Returns
        -------
        str
            The name for the heater costs.
        """
        return f"{self.central_or_decentral} resistive heater"

    @property
    def gas_boiler_costs_name(self) -> str:
        """
        Generates the name for the gas boiler costs based on the heat system.
        Used to retrieve data from `technology-data <https://github.com/PyPSA/technology-data>`.

        Returns
        -------
        str
            The name for the gas boiler costs.
        """
        return f"{self.central_or_decentral} gas boiler"

    @property
    def oil_boiler_costs_name(self) -> str:
        """
        Generates the name for the oil boiler costs based on the heat system.
        Used to retrieve data from `technology-data <https://github.com/PyPSA/technology-data>`.

        Returns
        -------
        str
            The name for the oil boiler costs.
        """
        return "decentral oil boiler"
