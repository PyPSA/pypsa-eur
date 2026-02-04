# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging
from enum import Enum

logger = logging.getLogger(__name__)


class HeatSourceType(Enum):
    """
    Categorization of heat sources by their fundamental characteristics.

    Attributes
    ----------
    INEXHAUSTIBLE : str
        Ambient sources (air, ground, sea water) that are always available
        and don't require resource tracking.
    SUPPLY_LIMITED : str
        Sources with spatial/temporal constraints requiring resource buses
        and generators (geothermal, river water).
    STORAGE : str
        Thermal storage discharge (PTES) with time-varying availability.
    PROCESS_WASTE : str
        Industrial process waste heat (PTX) that is a byproduct of other
        network components.
    """

    INEXHAUSTIBLE = "inexhaustible"
    SUPPLY_LIMITED = "supply_limited"
    STORAGE = "storage"
    PROCESS_WASTE = "process_waste"


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
    # PTX excess heat sources
    ELECTROLYSIS_waste = "electrolysis_waste"
    FISCHER_TROPSCH_waste = "fischer_tropsch_waste"
    SABATIER_waste = "sabatier_waste"
    HABER_BOSCH_waste = "haber_bosch_waste"
    METHANOLISATION_waste = "methanolisation_waste"
    FUEL_CELL_waste = "fuel_cell_waste"

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
    def source_type(self) -> HeatSourceType:
        """
        Get the category of this heat source.

        Returns
        -------
        HeatSourceType
            The category: INEXHAUSTIBLE, SUPPLY_LIMITED, STORAGE, or PROCESS_WASTE.
        """
        if self in [HeatSource.AIR, HeatSource.GROUND, HeatSource.SEA_WATER]:
            return HeatSourceType.INEXHAUSTIBLE
        elif self in [HeatSource.GEOTHERMAL, HeatSource.RIVER_WATER]:
            return HeatSourceType.SUPPLY_LIMITED
        elif self == HeatSource.PTES:
            return HeatSourceType.STORAGE
        else:
            return HeatSourceType.PROCESS_WASTE

    @property
    def temperature_from_config(self) -> bool:
        """
        Check if the heat source temperature is specified in config.

        Returns True if the temperature is a scalar value from config
        (heat_source_temperatures), False if it comes from a time-series file.

        Returns
        -------
        bool
            True for sources with config-defined temperatures (geothermal, PTX).
            False for sources with file-based time-series (river_water, ptes).
        """
        if self == HeatSource.RIVER_WATER:
            return False
        return self.source_type in [
            HeatSourceType.SUPPLY_LIMITED,
            HeatSourceType.PROCESS_WASTE,
        ]

    @property
    def requires_bus(self) -> bool:
        """
        Check whether the heat source requires a resource bus.

        Limited heat sources require a resource bus to track availability.
        Inexhaustible sources (air, ground, sea water) are always available
        and don't need resource tracking.

        Returns
        -------
        bool
            True if not INEXHAUSTIBLE, False otherwise.
        """
        return self.source_type != HeatSourceType.INEXHAUSTIBLE

    @property
    def requires_generator(self) -> bool:
        """
        Check if the heat source requires a generator component in the network.

        Generators represent the extraction capacity from geothermal wells or
        river water intakes, with associated capital costs and lifetimes.

        Returns
        -------
        bool
            True for SUPPLY_LIMITED sources only.
        """
        return self.source_type == HeatSourceType.SUPPLY_LIMITED

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
            True for STORAGE, PROCESS_WASTE, and GEOTHERMAL sources.
        """
        if self.source_type in [HeatSourceType.STORAGE, HeatSourceType.PROCESS_WASTE]:
            return True
        if self == HeatSource.GEOTHERMAL:
            return True
        else:
            return False

    @property
    def process_carrier(self) -> str | None:
        """
        Get the carrier name of the industrial process producing this waste heat.

        Returns
        -------
        str or None
            The carrier name (e.g., "Fischer-Tropsch"), or None if not a process waste source.
        """
        mapping = {
            HeatSource.FISCHER_TROPSCH_waste: "Fischer-Tropsch",
            HeatSource.SABATIER_waste: "Sabatier",
            HeatSource.HABER_BOSCH_waste: "Haber-Bosch",
            HeatSource.METHANOLISATION_waste: "methanolisation",
            HeatSource.ELECTROLYSIS_waste: "H2 Electrolysis",
            HeatSource.FUEL_CELL_waste: "H2 Fuel Cell",
        }
        return mapping.get(self)

    @property
    def process_output_bus_index(self) -> int | None:
        """
        Get which bus index (2, 3, or 4) the waste heat output uses on the process link.

        Returns
        -------
        int or None
            The bus index for efficiency/bus assignment, or None if not a process waste source.
        """
        mapping = {
            HeatSource.FISCHER_TROPSCH_waste: 3,
            HeatSource.SABATIER_waste: 3,
            HeatSource.HABER_BOSCH_waste: 3,
            HeatSource.METHANOLISATION_waste: 4,
            HeatSource.ELECTROLYSIS_waste: 2,
            HeatSource.FUEL_CELL_waste: 2,
        }
        return mapping.get(self)

    @property
    def waste_heat_option_key(self) -> str | None:
        """
        Get the config option key controlling this waste heat source utilisation.

        Returns
        -------
        str or None
            The option key (e.g., "use_fischer_tropsch_waste_heat"), or None if not applicable.
        """
        mapping = {
            HeatSource.FISCHER_TROPSCH_waste: "use_fischer_tropsch_waste_heat",
            HeatSource.SABATIER_waste: "use_methanation_waste_heat",
            HeatSource.HABER_BOSCH_waste: "use_haber_bosch_waste_heat",
            HeatSource.METHANOLISATION_waste: "use_methanolisation_waste_heat",
            HeatSource.ELECTROLYSIS_waste: "use_electrolysis_waste_heat",
            HeatSource.FUEL_CELL_waste: "use_fuel_cell_waste_heat",
        }
        return mapping.get(self)

    @property
    def technology_data_name(self) -> str | None:
        """
        Get the costs.csv technology name for efficiency-heat lookup.

        Some processes (Sabatier, Fuel Cell) use calculated efficiencies based on
        energy balance rather than a costs lookup.

        Returns
        -------
        str or None
            The technology name for costs lookup, or None if efficiency is calculated.
        """
        mapping = {
            HeatSource.FISCHER_TROPSCH_waste: "Fischer-Tropsch",
            HeatSource.HABER_BOSCH_waste: "Haber-Bosch",
            HeatSource.ELECTROLYSIS_waste: "electrolysis",
            HeatSource.HABER_BOSCH_waste: "Haber-Bosch",
            HeatSource.METHANOLISATION_waste: "methanolisation",
        }
        return mapping.get(self)

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

        For direct utilisation heat sources (geothermal), retrieves cost from technology-data.
        For other limited sources (like river_water without direct utilisation), returns 0.0.
        For inexhaustible sources (air, ground, sea water), this method shouldn't be called.

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

        For direct utilisation heat sources (geothermal), retrieves lifetime from technology-data.
        For other limited sources (like river_water without direct utilisation), returns infinity.
        For inexhaustible sources (air, ground, sea water), this method shouldn't be called.

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

    def get_waste_heat_efficiency(
        self, n, costs, nodes, fallback_ptx_heat_losses: float
    ):
        """
        Get the waste heat efficiency for this PTX process.

        For most processes, this looks up "efficiency-heat" from the costs data.
        For Sabatier and Fuel Cell, efficiency is calculated from energy balance
        as (1 - losses - main_efficiency).

        Parameters
        ----------
        n : pypsa.Network
            The network containing the process links.
        costs : pd.DataFrame
            DataFrame containing technology cost and efficiency parameters.
        nodes : pd.Index
            Node identifiers for looking up link efficiencies.

        Returns
        -------
        float or pd.Series
            The waste heat efficiency (heat output per unit of main input).

        Raises
        ------
        ValueError
            If called on a non-PROCESS_WASTE heat source.
        """
        if self.source_type != HeatSourceType.PROCESS_WASTE:
            raise ValueError(
                f"get_waste_heat_efficiency only applies to PROCESS_WASTE sources, not {self}"
            )

        if self == HeatSource.FISCHER_TROPSCH_waste:
            return costs.at[self.technology_data_name, "efficiency-heat"]
        elif self == HeatSource.ELECTROLYSIS_waste:
            return costs.at[self.technology_data_name, "efficiency-heat"]
        elif self == HeatSource.HABER_BOSCH_waste:
            return (
                costs.at[self.technology_data_name, "efficiency-heat"]
                / costs.at[self.technology_data_name, "electricity-input"]
            )
        elif self == HeatSource.METHANOLISATION_waste:
            return (
                costs.at[self.technology_data_name, "heat-output"]
                / costs.at[self.technology_data_name, "hydrogen-input"]
            )
        elif self == HeatSource.SABATIER_waste:
            return (
                1
                - fallback_ptx_heat_losses
                - n.links.loc[nodes + " Sabatier", "efficiency"]
            )
        elif self == HeatSource.FUEL_CELL_waste:
            return (
                1
                - fallback_ptx_heat_losses
                - n.links.loc[nodes + " H2 Fuel Cell", "efficiency"]
            )
        else:
            raise ValueError(f"No efficiency calculation defined for {self}")

    def get_heat_pump_input_bus(self, nodes, heat_system: str) -> str:
        """
        Get the secondary input bus for the heat pump link.

        The heat pump link has bus0 (electricity input), bus1 (heat output),
        and optionally bus2 (heat source input). This method determines bus2
        based on the heat source type:

        - Inexhaustible sources: No bus2 (empty string)
        - Limited sources: bus2 is the heat pump input bus (at return temperature if pre-heated, else at source temperature)

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
        if self.requires_bus:
            return nodes + f" {self.heat_pump_input_carrier(heat_system)}"
        else:
            return ""

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
            1 - 1/COP for limited sources (tracks heat drawn from source bus).

        Raises
        ------
        NotImplementedError
            If called for inexhaustible heat sources.

        See Also
        --------
        prepare_sector_network.add_heat : Creates heat pump links with this efficiency.
        """
        if self.requires_bus:
            return 1 - (1 / cop_heat_pump.clip(lower=0.001))
        else:
            return None

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

    def preheater_input_carrier(self, heat_system) -> str:
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
            Carrier name with '-pre-heater input' suffix in format '{heat_system} {source} pre-heater input'.
        """
        return f"{self.heat_carrier(heat_system)} pre-heater input"

    def heat_pump_input_carrier(self, heat_system) -> str:
        """
        Get the carrier name for fully-cooled heat from this source.

        Represents heat pump input, for
        final temperature lift by heat pump.

        Parameters
        ----------
        heat_system : HeatSystem or str
            The heat system (e.g., 'urban central').

        Returns
        -------
        str
            Carrier name with '-heat-pump input' suffix in format '{heat_system} {source} heat-pump input'.
        """
        return f"{self.heat_carrier(heat_system)} heat-pump input"

    def preheater_input_bus(self, nodes, heat_system) -> str:
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
            Bus name combining nodes with medium-temperature carrier in format 'nodes + {heat_system} {source} pre-heater input'.
        """
        return nodes + f" {self.preheater_input_carrier(heat_system)}"

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
