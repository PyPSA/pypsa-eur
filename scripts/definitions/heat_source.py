# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging
from enum import Enum

logger = logging.getLogger(__name__)


class HeatSource(Enum):
    """
    Enumeration representing different heat sources for heat pumps.

    Attributes
    ----------
    GEOTHERMAL : str
        Geothermal heat source.
    RIVER_WATER : str
        River water heat source.
    SEA_WATER : str
        Sea water heat source.
    AIR : str
        Air heat source.
    GROUND : str
        Ground heat source.
    PTES: str
        Pit Thermal Energy Storage as heat source.

    Methods
    -------
    __str__()
        Returns the string representation of the heat source.
    constant_temperature_celsius()
        Returns the constant temperature in Celsius for the heat source.
    is_limited()
        Returns whether the heat source is limited (vs. inexhaustible).
    requires_generator()
        Returns whether the heat source requires a generator.
    requires_preheater()
        Returns whether the heat source requires a preheater.
    supports_direct_utilisation()
        Returns whether the heat source supports direct heat utilisation.
    get_capital_cost(costs, overdim_factor, heat_system)
        Returns the capital cost for the heat source generator.
    get_lifetime(costs, heat_system)
        Returns the lifetime for the heat source generator.
    get_heat_pump_bus2(nodes, heat_system, heat_carrier)
        Returns the bus2 configuration for the heat pump link.
    get_heat_pump_efficiency2(cop_heat_pump)
        Returns the efficiency2 for the heat pump link.
    get_efficiency_pre_heater(efficiency_direct_utilisation)
        Returns the efficiency for the preheater link.
    get_efficiency_direct_utilisation(direct_heat_profile, nodes, n)
        Returns the efficiency for direct heat utilisation.
    """

    GEOTHERMAL = "geothermal"
    RIVER_WATER = "river_water"
    SEA_WATER = "sea_water"
    AIR = "air"
    GROUND = "ground"
    PTES = "ptes"

    def __init__(self, *args):
        super().__init__(*args)

    def __str__(self) -> str:
        """
        Returns the string representation of the heat source.

        Returns
        -------
        str
            The string representation of the heat source.
        """
        return self.value

    @property
    def has_constant_temperature(self) -> bool:
        """
        Returns the constant temperature in Celsius for the heat source.

        Returns
        -------
        bool
            True if the heat source has a constant temperature else False.
        """
        if self == HeatSource.GEOTHERMAL:
            return True
        else:
            return False

    @property
    def requires_bus(self) -> bool:
        """
        Returns whether the heat source is limited (vs. inexhaustible).

        Limited heat sources require a resource bus and have spatial/temporal constraints.
        Inexhaustible sources (air, ground) are always available.

        Returns
        -------
        bool
            True if the heat source is limited, False if inexhaustible.
        """
        if self in [
            HeatSource.GEOTHERMAL,
            HeatSource.RIVER_WATER,
        ]:
            return True
        else:
            return False

    @property
    def requires_generator(self) -> bool:
        """
        Returns whether the heat source requires a generator.

        Returns
        -------
        bool
            True if the heat source requires a generator, False otherwise.
        """
        if self in [HeatSource.GEOTHERMAL, HeatSource.RIVER_WATER]:
            return True
        else:
            return False

    @property
    def requires_preheater(self) -> bool:
        """
        Returns whether the heat source requires a pre-heater.

        Returns
        -------
        bool
            True if the heat source requires a pre-heater, False otherwise.
        """
        if self in [HeatSource.GEOTHERMAL, HeatSource.PTES]:
            return True
        else:
            return False

    def get_capital_cost(self, costs, overdim_factor: float, heat_system) -> float:
        """
        Returns the capital cost for the heat source generator.

        Parameters
        ----------
        costs : pd.DataFrame
            DataFrame containing cost information for different technologies.
        overdim_factor : float
            Factor to overdimension the heat generator.
        heat_system : HeatSystem
            The heat system for which to get the capital cost.

        Returns
        -------
        float
            The capital cost for the heat source generator.
        """
        # For direct utilisation heat sources, get cost from technology-data
        # For other limited sources (like river_water without direct utilisation), return 0.0
        # For inexhaustible sources, this method shouldn't be called
        if self in [HeatSource.GEOTHERMAL]:
            return (
                costs.at[
                    heat_system.heat_source_costs_name(self),
                    "capital_cost",
                ]
                * overdim_factor
            )
        else:
            return 0.0

    def get_lifetime(self, costs, heat_system) -> float:
        """
        Returns the lifetime for the heat source generator.

        Parameters
        ----------
        costs : pd.DataFrame
            DataFrame containing cost information for different technologies.
        heat_system : HeatSystem
            The heat system for which to get the lifetime.

        Returns
        -------
        float
            The lifetime for the heat source generator in years.
        """
        # For direct utilisation heat sources, get lifetime from technology-data
        # For other limited sources (like river_water without direct utilisation), return np.inf
        # For inexhaustible sources, this method shouldn't be called
        if self in [HeatSource.GEOTHERMAL]:
            return costs.at[heat_system.heat_source_costs_name(self), "lifetime"]
        else:
            return float("inf")

    def get_heat_pump_bus2(self, nodes, heat_carrier: str) -> str:
        """
        Returns the bus2 configuration for the heat pump link.

        Parameters
        ----------
        nodes : pd.Index or list
            The nodes for which to generate the bus name.
        requires_preheater : bool
            Whether the heat source requires a preheater.
        heat_carrier : str
            The heat carrier name.

        Returns
        -------
        str
            The bus2 name for the heat pump, or empty string if not applicable.
        """
        if self in [HeatSource.AIR, HeatSource.GROUND, HeatSource.SEA_WATER]:
            # Inexhaustible sources (air, ground, sea-water) don't have a bus2
            return ""
        elif self.requires_preheater:
            # Sources with preheater use pre-chilled bus
            return nodes + f" {heat_carrier} pre-chilled"
        else:
            # Limited sources without preheater use the heat carrier bus directly
            return nodes + f" {heat_carrier}"

    def get_heat_pump_efficiency2(self, cop_heat_pump) -> float:
        """
        Returns the efficiency2 for the heat pump link.

        Parameters
        ----------
        cop_heat_pump : float or pd.Series
            The coefficient of performance of the heat pump.

        Returns
        -------
        float or pd.Series
            The efficiency2 value for the heat pump.
        """
        if self in [HeatSource.AIR, HeatSource.GROUND, HeatSource.SEA_WATER]:
            # Inexhaustible sources (air, ground, sea-water) don't have a bus2
            # Inexhaustible sources (air, ground, sea_water) have efficiency2 = 1
            # (no resource consumption from dummy bus)
            return 1.0
        else:
            # Limited heat sources consume heat from the resource bus (either pre-chilled or direct)
            # efficiency2 represents the fraction of heat drawn from the source
            # This is 1 - (1/COP), representing (COP-1)/COP
            return 1 - (1 / cop_heat_pump.clip(lower=0.001))

    def requires_heat_pump(self, ptes_discharge_resistive_boosting: bool) -> bool:
        """
        Returns whether the heat source requires a heat pump.

        Returns
        -------
        bool
            True if the heat source requires a heat pump, False otherwise.
        """

        if self == HeatSource.PTES and ptes_discharge_resistive_boosting:
            logging.info(
                "PTES configured with resistive boosting during discharge; "
                "heat pump not built for PTES."
            )
            return False
        else:
            return True
