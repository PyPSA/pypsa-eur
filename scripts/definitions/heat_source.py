# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging
from enum import Enum

logger = logging.getLogger(__name__)


class HeatSource(Enum):
    """
    Enumeration representing different heat sources for heat pumps and direct utilisation.

    Heat sources are categorized by their characteristics:

    **Inexhaustible sources** (AIR, GROUND, SEA_WATER):
        Always available, no resource bus needed, heat pump draws from ambient.

    **Limited sources requiring a bus** (GEOTHERMAL, RIVER_WATER, PTES):
        Have spatial/temporal constraints, require resource tracking via buses.
        May support direct utilisation or preheating depending on temperature.

    **Sources with boosting** (PTES, geothermal, river water):
        When source temperature is below forward temperature, a heat pump
        provides supplemental boost energy routed through the hp_output_bus.

    Attributes
    ----------
    GEOTHERMAL : str
        Geothermal heat source with constant temperature.
    RIVER_WATER : str
        River water heat source with time-varying temperature.
    SEA_WATER : str
        Sea water heat source (treated as inexhaustible).
    AIR : str
        Ambient air heat source (inexhaustible).
    GROUND : str
        Ground/soil heat source (inexhaustible).
    PTES : str
        Pit Thermal Energy Storage discharge as heat source.

    See Also
    --------
    HeatSystem : Defines heat system types (urban central, urban decentral, rural).
    build_heat_source_utilisation_profiles : Calculates utilisation profiles for heat sources.
    build_cop_profiles : Calculates COP profiles for heat pumps using these sources.
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
        Check if the heat source has a constant (time-invariant) temperature.

        Constant-temperature sources (e.g., geothermal) have their temperature
        specified in config rather than loaded from time series files.

        Returns
        -------
        bool
            True for geothermal, False for all other sources.
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
            HeatSource.PTES,
        ]:
            return True
        else:
            return False

    @property
    def requires_generator(self) -> bool:
        """
        Check if the heat source requires a generator component in the network.

        Generators represent the extraction capacity from geothermal wells or
        river water intakes, with associated capital costs and lifetimes.

        Returns
        -------
        bool
            True for geothermal and river_water, False otherwise.
        """
        if self in [HeatSource.GEOTHERMAL, HeatSource.RIVER_WATER]:
            return True
        else:
            return False

    def requires_heat_pump(self, ptes_discharge_resistive_boosting: bool) -> bool:
        """
        Check if a heat pump should be built for this heat source.

        Most heat sources require a heat pump to lift temperature to the
        forward temperature. PTES is special: it can use either a heat pump
        or resistive boosting for temperature lift during discharge.

        Parameters
        ----------
        ptes_discharge_resistive_boosting : bool
            Whether PTES uses resistive heaters instead of heat pumps.

        Returns
        -------
        bool
            False for PTES with resistive boosting, True otherwise.
        """
        if self == HeatSource.PTES and ptes_discharge_resistive_boosting:
            logging.info(
                "PTES configured with resistive boosting during discharge; "
                "heat pump not built for PTES."
            )
            return False
        else:
            return True

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

        Notes
        -----
        - For direct utilisation heat sources (geothermal), gets cost from technology-data.
        - For other limited sources (like river_water), returns 0.0.
        - For inexhaustible sources, this method shouldn't be called.
        """
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

        Notes
        -----
        - For direct utilisation heat sources (geothermal), gets lifetime from technology-data.
        - For other limited sources (like river_water), returns infinity.
        - For inexhaustible sources, this method shouldn't be called.
        """
        if self in [HeatSource.GEOTHERMAL]:
            return costs.at[heat_system.heat_source_costs_name(self), "lifetime"]
        else:
            return float("inf")

    def heat_carrier(self, heat_system) -> str:
        """
        Get the carrier name for heat from this source.

        Parameters
        ----------
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Carrier name in format '{heat_system} {source} heat',
            e.g., 'urban central ptes heat'.
        """
        return f"{heat_system} {self} heat"

    def hp_input_carrier(self, heat_system) -> str:
        """
        Get the carrier name for heat produced by the boosting heat pump.

        This bus receives HP output and is consumed by the utilisation link
        at rate ``boosting_ratio`` per unit of source heat.

        Parameters
        ----------
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Carrier name in format '{heat_system} {source} heat hp output'.
        """
        return f"{self.heat_carrier(heat_system)} heat pump input"

    def hp_input_bus(self, nodes, heat_system) -> str:
        """
        Get the bus name for heat pump output at the given nodes.

        Parameters
        ----------
        nodes : pd.Index or str
            Node identifier(s).
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Bus name in format 'nodes + {heat_system} {source} heat hp output'.
        """
        if self.requires_bus:
            return nodes + f" {self.hp_input_carrier(heat_system)}"
        else:
            return nodes + f" {heat_system} heat"

    def resource_bus(self, nodes, heat_system) -> str:
        """
        Get the primary resource bus for heat from this source, e.g. `DE0 0 urban central geothermal heat` or `DE0 0 urban central ptes heat`.

        This is where heat enters the system from generators or storage
        discharge, at the source's native temperature.

        Parameters
        ----------
        nodes : pd.Index or str
            Node identifier(s).
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Bus name combining nodes with heat carrier in format 'nodes + {heat_system} {source} heat'.
        """
        return nodes + f" {self.heat_carrier(heat_system)}"
