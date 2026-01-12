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

    **Sources with preheater** (PTES):
        When source temperature is between return and forward temperatures,
        can preheat return flow before heat pump provides final lift.

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

    @property
    def requires_preheater(self) -> bool:
        """
        Check if the heat source uses preheating when below forward temperature.

        Preheating allows intermediate-temperature sources to warm the return
        flow before a heat pump provides the final temperature lift, improving
        overall system efficiency.

        Returns
        -------
        bool
            True for PTES and GEOTHERMAL, False otherwise.
        """
        if self in [HeatSource.PTES, HeatSource.GEOTHERMAL]:
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

    def get_heat_pump_bus2(self, nodes, heat_system: str) -> str:
        """
        Get the secondary input bus for the heat pump link.

        The heat pump link has bus0 (electricity input), bus1 (heat output),
        and optionally bus2 (heat source input). This method determines bus2
        based on the heat source type:

        - Inexhaustible sources: No bus2 (empty string)
        - Sources with preheater: Return-temperature bus (post-preheat)
        - Other limited sources: Resource bus directly

        Parameters
        ----------
        nodes : pd.Index or list
            The nodes for which to generate the bus name.
        heat_system : str
            The heat system identifier (e.g., 'urban central').

        Returns
        -------
        str
            The bus2 name for the heat pump, or empty string if not applicable.
        """
        if self in [HeatSource.AIR, HeatSource.GROUND, HeatSource.SEA_WATER]:
            return ""
        elif self.requires_preheater:
            return self.return_temperature_bus(nodes, heat_system)
        else:
            return self.resource_bus(nodes, heat_system)

    def get_heat_pump_efficiency2(self, cop_heat_pump) -> float:
        """
        Get the efficiency2 (heat source consumption) for the heat pump link.

        The heat pump link uses an inverted bus configuration where bus0 is the
        heat output (not input) to attribute capital costs to the heat bus. This
        means the link operates with negative flow (p_max_pu=0, p_min_pu < 0).

        - Inexhaustible sources (air, ground, sea_water): Returns 1.0 since no
          resource tracking is needed.
        - Limited sources (geothermal, river_water, ptes): Returns 1 - 1/COP to
          track heat drawn from the source bus.

        For a standard heat pump energy balance:

            Q_output = Q_source + W_electricity
            COP = Q_output / W_electricity

        The fraction of output heat from the source is:

            Q_source / Q_output = (COP - 1) / COP

        However, since bus0 is the output in our inverted configuration, PyPSA
        interprets efficiency2 as: input_bus2 = p_bus0 * efficiency2

        With negative p_bus0 (heat flowing out), we need efficiency2 to give the
        correct heat drawn from the source. The formula simplifies to:

            efficiency2 = 1 - 1/COP

        This ensures that for each unit of heat output, the link draws the correct
        amount from the heat source bus.

        Parameters
        ----------
        cop_heat_pump : float or pd.Series
            The coefficient of performance of the heat pump.

        Returns
        -------
        float or pd.Series
            1.0 for inexhaustible sources (no resource tracking needed),
            1 - 1/COP for limited sources (tracks heat drawn from source bus).

        See Also
        --------
        prepare_sector_network.add_heat : Creates heat pump links with this efficiency.
        """
        if self in [HeatSource.AIR, HeatSource.GROUND, HeatSource.SEA_WATER]:
            return 1.0
        else:
            return 1 - (1 / cop_heat_pump.clip(lower=0.001))

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

    def medium_temperature_carrier(self, heat_system) -> str:
        """
        Get the carrier name for partially-cooled heat from this source.

        Used in cascading temperature utilisation: heat that has been used
        for direct supply but still has usable thermal energy.

        Parameters
        ----------
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Carrier name with '-medium-temperature' suffix in format '{heat_system} {source} medium-temperature'.
        """
        return f"{self.heat_carrier(heat_system)} medium-temperature"

    def return_temperature_carrier(self, heat_system) -> str:
        """
        Get the carrier name for fully-cooled heat from this source.

        Represents heat at return temperature after preheating, ready for
        final temperature lift by heat pump.

        Parameters
        ----------
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Carrier name with '-return-temperature' suffix in format '{heat_system} {source} return-temperature'.
        """
        return f"{self.heat_carrier(heat_system)} return-temperature"

    def medium_temperature_bus(self, nodes, heat_system) -> str:
        """
        Get bus name for partially-cooled heat at the given nodes.

        Parameters
        ----------
        nodes : pd.Index or str
            Node identifier(s).
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Bus name combining nodes with medium-temperature carrier in format 'nodes + {heat_system} {source} medium-temperature'.
        """
        return nodes + f" {self.medium_temperature_carrier(heat_system)}"

    def return_temperature_bus(self, nodes, heat_system) -> str:
        """
        Get bus name for fully-cooled heat at the given nodes.

        Parameters
        ----------
        nodes : pd.Index or str
            Node identifier(s).
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Bus name combining nodes with return-temperature carrier in format 'nodes + {heat_system} {source} return-temperature'.
        """
        return nodes + f" {self.return_temperature_carrier(heat_system)}"

    def resource_bus(self, nodes, heat_system) -> str:
        """
        Get the primary resource bus for heat from this source.

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
