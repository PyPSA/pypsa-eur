# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Sector configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#sector
"""

from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from scripts.lib.validation.config._base import ConfigModel


class _DistrictHeatingConfig(ConfigModel):
    """Configuration for `sector.district_heating` settings."""

    potential: float | dict[str, float] = Field(
        0.6,
        description="Maximum fraction of urban demand which can be supplied by district heating. If given as dictionary, specify one value per country modeled or provide a default value with key `default` to fill values for all unspecified countries.",
    )
    progress: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.0,
            2025: 0.1,
            2030: 0.25,
            2035: 0.4,
            2040: 0.55,
            2045: 0.75,
            2050: 1.0,
        },
        description="Increase of today's district heating demand to potential maximum district heating share. Progress = 0 means today's district heating share. Progress = 1 means maximum fraction of urban demand is supplied by district heating.",
    )
    district_heating_loss: float = Field(
        0.15,
        description="Share increase in district heat demand in urban central due to heat losses.",
    )
    supply_temperature_approximation: dict[str, Any] = Field(
        default_factory=lambda: {
            "max_forward_temperature_baseyear": {
                "FR": 110,
                "DK": 75,
                "DE": 109,
                "CZ": 130,
                "FI": 115,
                "PL": 130,
                "SE": 102,
                "IT": 90,
            },
            "min_forward_temperature_baseyear": {"DE": 82},
            "return_temperature_baseyear": {"DE": 58},
            "lower_threshold_ambient_temperature": 0,
            "upper_threshold_ambient_temperature": 10,
            "rolling_window_ambient_temperature": 72,
            "relative_annual_temperature_reduction": 0.01,
        },
        description="Supply temperature approximation settings.",
    )
    ptes: dict[str, Any] = Field(
        default_factory=lambda: {
            "dynamic_capacity": False,
            "supplemental_heating": {"enable": False, "booster_heat_pump": False},
            "max_top_temperature": 90,
            "min_bottom_temperature": 35,
        },
        description="Pit thermal energy storage settings.",
    )
    ates: dict[str, Any] = Field(
        default_factory=lambda: {
            "enable": False,
            "suitable_aquifer_types": ["Highly productive porous aquifers"],
            "aquifer_volumetric_heat_capacity": 2600,
            "fraction_of_aquifer_area_available": 0.2,
            "effective_screen_length": 20,
            "capex_as_fraction_of_geothermal_heat_source": 0.75,
            "recovery_factor": 0.6,
            "marginal_cost_charger": 0.035,
            "ignore_missing_regions": False,
        },
        description="Aquifer thermal energy storage settings.",
    )
    heat_source_cooling: float = Field(
        6, description="Cooling of heat source for heat pumps."
    )
    heat_pump_cop_approximation: dict[str, Any] = Field(
        default_factory=lambda: {
            "refrigerant": "ammonia",
            "heat_exchanger_pinch_point_temperature_difference": 5,
            "isentropic_compressor_efficiency": 0.8,
            "heat_loss": 0.0,
            "min_delta_t_lift": 10,
        },
        description="Heat pump COP approximation settings.",
    )
    limited_heat_sources: dict[str, Any] = Field(
        default_factory=lambda: {
            "geothermal": {
                "constant_temperature_celsius": 65,
                "ignore_missing_regions": False,
            },
            "river_water": {"constant_temperature_celsius": False},
        },
        description="Dictionary with names of limited heat sources (not air). Must be `river_water` / `geothermal` or another heat source in `Manz et al. 2024 <https://www.sciencedirect.com/science/article/pii/S0960148124001769>`_.",
    )
    direct_utilisation_heat_sources: list[str] = Field(
        default_factory=lambda: ["geothermal"],
        description="List of heat sources for direct heat utilisation in district heating. Must be in the keys of `heat_utilisation_potentials` (e.g. `geothermal`).",
    )
    temperature_limited_stores: list[str] = Field(
        default_factory=lambda: ["ptes"],
        description="List of names for stores used as limited heat sources.",
    )
    dh_areas: dict[str, Any] = Field(
        default_factory=lambda: {"buffer": 1000, "handle_missing_countries": "fill"},
        description="District heating areas settings.",
    )


class _ResidentialHeatDsmConfig(BaseModel):
    """Configuration for `sector.residential_heat.dsm` settings."""

    enable: bool = Field(
        False,
        description="Enable residential heat demand-side management that allows heating systems to provide flexibility by shifting demand within configurable time periods. Models building thermal mass as energy storage.",
    )
    direction: list[str] = Field(
        default_factory=lambda: ["overheat", "undercool"],
        description="'overheat-undercool' means both pre-heating and delayed heating are allowed. 'overheat' allows only pre-heating where buildings are heated up above target temperature and then allowed to cool down, while 'undercool' allows only delayed heating where buildings can cool below target temperature and then be heated up again.",
    )
    restriction_value: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.06,
            2025: 0.16,
            2030: 0.27,
            2035: 0.36,
            2040: 0.38,
            2045: 0.39,
            2050: 0.4,
        },
        description="Maximum state of charge (as fraction) for heat flexibility storage representing available thermal buffer capacity in buildings. Set to 0 for no flexibility or to 1.0 to assume that the entire heating demand can contribute to flexibility.",
    )
    restriction_time: list[int] = Field(
        default_factory=lambda: [10, 22],
        description="Checkpoint hours (0-23) at which heat flexibility storage must return to baseline state of charge, i.e. the residence surplus or missing heat be balanced. Time is the local time for each country and bus. Default: [10, 22] creates 12-hour periods with checkpoints at 10am and 10pm.",
    )


class _ResidentialHeatConfig(BaseModel):
    """Configuration for `sector.residential_heat` settings."""

    dsm: _ResidentialHeatDsmConfig = Field(
        default_factory=_ResidentialHeatDsmConfig,
        description="Configuration options for residential heat demand-side management (DSM). See `smartEn DSM study <https://smarten.eu/wp-content/uploads/2022/09/SmartEn-DSF-benefits-2030-Report_DIGITAL.pdf>`_ (Appendix A) for methodology.",
    )


class _RetrofittingConfig(BaseModel):
    """Configuration for `sector.retrofitting` settings."""

    retro_endogen: bool = Field(
        False,
        description="Add retrofitting as an endogenous system which co-optimise space heat savings.",
    )
    cost_factor: float = Field(1.0, description="Weight costs for building renovation.")
    interest_rate: float = Field(
        0.04, description="The interest rate for investment in building components."
    )
    annualise_cost: bool = Field(
        True, description="Annualise the investment costs of retrofitting."
    )
    tax_weighting: bool = Field(
        False,
        description="Weight the costs of retrofitting depending on taxes in countries.",
    )
    construction_index: bool = Field(
        True,
        description="Weight the costs of retrofitting depending on labour/material costs per country.",
    )


class _CHPConfig(BaseModel):
    """Configuration for `sector.chp` settings."""

    enable: bool = Field(
        True, description="Add option for using Combined Heat and Power (CHP)."
    )
    fuel: list[str] = Field(
        default_factory=lambda: ["solid biomass", "gas"],
        description='Possible options are all fuels which have an existing bus and their CO2 intensity is given in the technology data. Currently possible are "gas", "oil", "methanol", "lignite", "coal" as well as "solid biomass". For all fuels except solid biomass, the techno-economic data from gas CHP is used. For the special case of solid biomass fuel, both CHP plants with and without carbon capture are added.',
    )
    micro_chp: bool = Field(
        False,
        description="Add option for using gas-fired Combined Heat and Power (CHP) for decentral areas.",
    )


class _MethanolConfig(BaseModel):
    """Configuration for `sector.methanol` settings."""

    regional_methanol_demand: bool = Field(
        False,
        description="Spatially resolve methanol demand. Set to true if regional CO2 constraints needed.",
    )
    methanol_reforming: bool = Field(False, description="Add methanol reforming.")
    methanol_reforming_cc: bool = Field(
        False, description="Add methanol reforming with carbon capture."
    )
    methanol_to_kerosene: bool = Field(False, description="Add methanol to kerosene.")
    methanol_to_power: dict[str, bool] = Field(
        default_factory=lambda: {
            "ccgt": False,
            "ccgt_cc": False,
            "ocgt": True,
            "allam": False,
        },
        description="Add different methanol to power technologies.",
    )
    biomass_to_methanol: bool = Field(True, description="Add biomass to methanol.")
    biomass_to_methanol_cc: bool = Field(
        False, description="Add biomass to methanol with carbon capture."
    )


class _TransmissionEfficiencyConfig(BaseModel):
    """Configuration for `sector.transmission_efficiency` settings."""

    enable: list[str] = Field(
        default_factory=lambda: [
            "DC",
            "H2 pipeline",
            "gas pipeline",
            "electricity distribution grid",
        ],
        description="Switch to select the carriers for which transmission efficiency is to be added. Carriers not listed assume lossless transmission.",
    )
    DC: dict[str, float] = Field(
        default_factory=lambda: {
            "efficiency_static": 0.98,
            "efficiency_per_1000km": 0.977,
        },
        description="DC transmission efficiency.",
    )
    H2_pipeline: dict[str, float] = Field(
        default_factory=lambda: {
            "efficiency_per_1000km": 1,
            "compression_per_1000km": 0.018,
        },
        alias="H2 pipeline",
        description="H2 pipeline transmission efficiency.",
    )
    gas_pipeline: dict[str, float] = Field(
        default_factory=lambda: {
            "efficiency_per_1000km": 1,
            "compression_per_1000km": 0.01,
        },
        alias="gas pipeline",
        description="Gas pipeline transmission efficiency.",
    )
    electricity_distribution_grid: dict[str, float] = Field(
        default_factory=lambda: {"efficiency_static": 0.97},
        alias="electricity distribution grid",
        description="Electricity distribution grid efficiency.",
    )

    model_config = ConfigDict(populate_by_name=True)


class _LimitMaxGrowthConfig(BaseModel):
    """Configuration for `sector.limit_max_growth` settings."""

    enable: bool = Field(
        False, description="Add option to limit the maximum growth of a carrier."
    )
    factor: float = Field(
        1.3,
        description="The maximum growth factor of a carrier (e.g. 1.3 allows  30% larger than max historic growth).",
    )
    max_growth: dict[str, float] = Field(
        default_factory=lambda: {
            "onwind": 16,
            "solar": 28,
            "offwind-ac": 35,
            "offwind-dc": 35,
        },
        description="The historic maximum growth of a carrier.",
    )
    max_relative_growth: dict[str, float] = Field(
        default_factory=lambda: {
            "onwind": 3,
            "solar": 3,
            "offwind-ac": 3,
            "offwind-dc": 3,
        },
        description="The historic maximum relative growth of a carrier.",
    )


class _EnhancedGeothermalConfig(BaseModel):
    """Configuration for `sector.enhanced_geothermal` settings."""

    enable: bool = Field(
        False, description="Add option to include Enhanced Geothermal Systems."
    )
    flexible: bool = Field(
        True, description="Add option for flexible operation (see Ricks et al. 2024)."
    )
    max_hours: int = Field(
        240,
        description="The maximum hours the reservoir can be charged under flexible operation.",
    )
    max_boost: float = Field(
        0.25, description="The maximum boost in power output under flexible operation."
    )
    var_cf: bool = Field(
        True,
        description="Add option for variable capacity factor (see Ricks et al. 2024).",
    )
    sustainability_factor: float = Field(
        0.0025,
        description="Share of sourced heat that is replenished by the earth's core (see details in `build_egs_potentials.py <https://github.com/PyPSA/pypsa-eur-sec/blob/master/scripts/build_egs_potentials.py>`_).",
    )


class _SolidBiomassImportConfig(BaseModel):
    """Configuration for `sector.solid_biomass_import` settings."""

    enable: bool = Field(
        False, description="Add option to include solid biomass imports."
    )
    price: float = Field(
        54, description="Price for importing solid biomass (currency/MWh)."
    )
    max_amount: float = Field(
        1390, description="Maximum solid biomass import potential (TWh)."
    )
    upstream_emissions_factor: float = Field(
        0.1, description="Upstream emissions of solid biomass imports."
    )


class _ImportsConfig(BaseModel):
    """Configuration for `sector.imports` settings."""

    enable: bool = Field(
        False, description="Add option to include renewable energy imports."
    )
    limit: float = Field(
        float("inf"), description="Maximum allowed renewable energy imports (TWh)."
    )
    limit_sense: str = Field("<=", description="Sense of the limit.")
    price: dict[str, float] = Field(
        default_factory=lambda: {
            "H2": 74,
            "NH3": 97,
            "methanol": 121,
            "gas": 122,
            "oil": 125,
        },
        description="Price for importing renewable energy of carrier.",
    )


class SectorConfig(BaseModel):
    """Configuration for `sector` settings."""

    transport: bool = Field(True, description="Flag to include transport sector.")
    heating: bool = Field(True, description="Flag to include heating sector.")
    biomass: bool = Field(True, description="Flag to include biomass sector.")
    industry: bool = Field(True, description="Flag to include industry sector.")
    shipping: bool = Field(True, description="Flag to include shipping sector.")
    aviation: bool = Field(True, description="Flag to include aviation sector.")
    agriculture: bool = Field(True, description="Flag to include agriculture sector.")
    fossil_fuels: bool = Field(
        True, description="Flag to include imports of fossil fuels."
    )

    district_heating: _DistrictHeatingConfig = Field(
        default_factory=_DistrictHeatingConfig,
        description="District heating configuration.",
    )

    heat_pump_sources: dict[str, list[str]] = Field(
        default_factory=lambda: {
            "urban central": ["air"],
            "urban decentral": ["air"],
            "rural": ["air", "ground"],
        },
        description="Heat pump sources by area.",
    )

    residential_heat: _ResidentialHeatConfig = Field(
        default_factory=_ResidentialHeatConfig,
        description="Residential heat configuration.",
    )

    cluster_heat_buses: bool = Field(
        True,
        description="Cluster residential and service heat buses in `prepare_sector_network.py <https://github.com/PyPSA/pypsa-eur-sec/blob/master/scripts/prepare_sector_network.py>`_ to one to save memory.",
    )
    heat_demand_cutout: str = Field("default", description="Heat demand cutout.")

    # Transport settings
    bev_dsm_restriction_value: float = Field(
        0.8,
        description="Adds a lower state of charge (SOC) limit for battery electric vehicles (BEV) to manage its own energy demand (DSM). Located in `build_transport_demand.py <https://github.com/PyPSA/pypsa-eur-sec/blob/master/scripts/build_transport_demand.py>`_. Set to 0 for no restriction on BEV DSM.",
    )
    bev_dsm_restriction_time: float = Field(
        7, description="Time at which SOC of BEV has to be dsm_restriction_value."
    )
    transport_heating_deadband_upper: float = Field(
        20.0,
        description="The maximum temperature in the vehicle. At higher temperatures, the energy required for cooling in the vehicle increases.",
    )
    transport_heating_deadband_lower: float = Field(
        15.0,
        description="The minimum temperature in the vehicle. At lower temperatures, the energy required for heating in the vehicle increases.",
    )
    ICE_lower_degree_factor: float = Field(
        0.375,
        description="Share increase in energy demand in internal combustion engine (ICE) for each degree difference between the cold environment and the minimum temperature.",
    )
    ICE_upper_degree_factor: float = Field(
        1.6,
        description="Share increase in energy demand in internal combustion engine (ICE) for each degree difference between the hot environment and the maximum temperature.",
    )
    EV_lower_degree_factor: float = Field(
        0.98,
        description="Share increase in energy demand in electric vehicles (EV) for each degree difference between the cold environment and the minimum temperature.",
    )
    EV_upper_degree_factor: float = Field(
        0.63,
        description="Share increase in energy demand in electric vehicles (EV) for each degree difference between the hot environment and the maximum temperature.",
    )
    bev_dsm: bool = Field(
        True,
        description="Add the option for battery electric vehicles (BEV) to participate in demand-side management (DSM).",
    )
    bev_dsm_availability: float = Field(
        0.5,
        description="The share for battery electric vehicles (BEV) that are able to do demand side management (DSM).",
    )
    bev_energy: float = Field(
        0.05, description="The average size of battery electric vehicles (BEV) in MWh."
    )
    bev_charge_efficiency: float = Field(
        0.9,
        description="Battery electric vehicles (BEV) charge and discharge efficiency.",
    )
    bev_charge_rate: float = Field(
        0.011,
        description="The power consumption for one electric vehicle (EV) in MWh. Value derived from 3-phase charger with 11 kW.",
    )
    bev_avail_max: float = Field(
        0.95,
        description="The maximum share plugged-in availability for passenger electric vehicles.",
    )
    bev_avail_mean: float = Field(
        0.8,
        description="The average share plugged-in availability for passenger electric vehicles.",
    )
    v2g: bool = Field(
        True,
        description="Allows feed-in to grid from EV battery. This is only enabled if BEV demand-side management is enabled, and the share of vehicles participating is V2G is given by `bev_dsm_availability`.",
    )

    land_transport_fuel_cell_share: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0,
            2025: 0,
            2030: 0,
            2035: 0,
            2040: 0,
            2045: 0,
            2050: 0,
        },
        description="The share of vehicles that uses fuel cells in a given year.",
    )
    land_transport_electric_share: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0,
            2025: 0.05,
            2030: 0.2,
            2035: 0.45,
            2040: 0.7,
            2045: 0.85,
            2050: 1,
        },
        description="The share of vehicles that uses electric vehicles (EV) in a given year.",
    )
    land_transport_ice_share: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 1,
            2025: 0.95,
            2030: 0.8,
            2035: 0.55,
            2040: 0.3,
            2045: 0.15,
            2050: 0,
        },
        description="The share of vehicles that uses internal combustion engines (ICE) in a given year. What is not EV or FCEV is oil-fuelled ICE.",
    )

    transport_electric_efficiency: float = Field(
        53.19,
        description="The conversion efficiencies of electric vehicles in transport.",
    )
    transport_fuel_cell_efficiency: float = Field(
        30.003, description="The H2 conversion efficiencies of fuel cells in transport."
    )
    transport_ice_efficiency: float = Field(
        16.0712,
        description="The oil conversion efficiencies of internal combustion engine (ICE) in transport.",
    )

    agriculture_machinery_electric_share: float = Field(
        0.5, description="The share for agricultural machinery that uses electricity."
    )
    agriculture_machinery_oil_share: float = Field(
        0.5, description="The share for agricultural machinery that uses oil."
    )
    agriculture_machinery_fuel_efficiency: float = Field(
        0.7,
        description="The efficiency of electric-powered machinery in the conversion of electricity to meet agricultural needs.",
    )
    agriculture_machinery_electric_efficiency: float = Field(
        0.3,
        description="The efficiency of oil-powered machinery in the conversion of oil to meet agricultural needs.",
    )

    shipping_hydrogen_liquefaction: bool = Field(
        False,
        description="Whether to include liquefaction costs for hydrogen demand in shipping.",
    )
    shipping_hydrogen_share: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0,
            2025: 0,
            2030: 0,
            2035: 0,
            2040: 0,
            2045: 0,
            2050: 0,
        },
        description="The share of ships powered by hydrogen in a given year.",
    )
    shipping_methanol_share: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0,
            2025: 0,
            2030: 0.15,
            2035: 0.35,
            2040: 0.55,
            2045: 0.8,
            2050: 1,
        },
        description="The share of ships powered by methanol in a given year.",
    )
    shipping_oil_share: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 1,
            2025: 1,
            2030: 0.85,
            2035: 0.65,
            2040: 0.45,
            2045: 0.2,
            2050: 0,
        },
        description="The share of ships powered by oil in a given year.",
    )
    shipping_methanol_efficiency: float = Field(
        0.46,
        description="The efficiency of methanol-powered ships in the conversion of methanol to meet shipping needs (propulsion). The efficiency increase from oil can be 10-15% higher according to the `IEA <https://www.iea-amf.org/app/webroot/files/file/Annex%20Reports/AMF_Annex_56.pdf>`_.",
    )
    shipping_oil_efficiency: float = Field(
        0.40,
        description="The efficiency of oil-powered ships in the conversion of oil to meet shipping needs (propulsion). Base value derived from 2011.",
    )

    aviation_demand_factor: float = Field(
        1.0,
        description="The proportion of demand for aviation compared to today's consumption.",
    )
    HVC_demand_factor: float = Field(
        1.0,
        description="The proportion of demand for high-value chemicals compared to today's consumption.",
    )

    time_dep_hp_cop: bool = Field(
        True,
        description="Consider the time dependent coefficient of performance (COP) of the heat pump.",
    )
    heat_pump_sink_T_individual_heating: float = Field(
        55.0,
        description="The temperature heat sink used in heat pumps based on DTU / large area radiators. The value is conservatively high to cover hot water and space heating in poorly-insulated buildings.",
    )

    reduce_space_heat_exogenously: bool = Field(
        True,
        description="Influence on space heating demand by a certain factor (applied before losses in district heating).",
    )
    reduce_space_heat_exogenously_factor: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0.10,
            2025: 0.09,
            2030: 0.09,
            2035: 0.11,
            2040: 0.16,
            2045: 0.21,
            2050: 0.29,
        },
        description="A positive factor can mean renovation or demolition of a building. If the factor is negative, it can mean an increase in floor area, increased thermal comfort, population growth. The default factors are determined by the `Eurocalc Homes and buildings decarbonization scenario <http://tool.european-calculator.eu/app/buildings/building-types-area/?levers=1ddd4444421213bdbbbddd44444ffffff11f411111221111211l212221>`_.",
    )

    retrofitting: _RetrofittingConfig = Field(
        default_factory=_RetrofittingConfig, description="Retrofitting configuration."
    )

    tes: bool = Field(
        True,
        description="Add option for storing thermal energy in large water pits associated with district heating systems and individual thermal energy storage (TES).",
    )
    boilers: bool = Field(
        True, description="Add option for transforming gas into heat using gas boilers."
    )
    resistive_heaters: bool = Field(
        True,
        description="Add option for transforming electricity into heat using resistive heaters (independently from gas boilers).",
    )
    oil_boilers: bool = Field(
        False, description="Add option for transforming oil into heat using boilers."
    )
    biomass_boiler: bool = Field(
        True, description="Add option for transforming biomass into heat using boilers."
    )
    overdimension_heat_generators: dict[str, float] = Field(
        default_factory=lambda: {"decentral": 1.1, "central": 1.0},
        description="Add option for overdimensioning heating systems by a certain factor. This allows them to cover heat demand peaks e.g. 10% higher than those in the data with a setting of 1.1.",
    )

    chp: _CHPConfig = Field(
        default_factory=_CHPConfig, description="CHP configuration."
    )
    solar_thermal: bool = Field(
        True, description="Add option for using solar thermal to generate heat."
    )
    solar_cf_correction: float = Field(
        0.788457,
        description="The correction factor for the value provided by the solar thermal profile calculations.",
    )

    methanation: bool = Field(
        True,
        description="Add option for transforming hydrogen and CO2 into methane using methanation.",
    )
    coal_cc: bool = Field(
        False, description="Add option for coal CHPs with carbon capture."
    )
    dac: bool = Field(True, description="Add option for Direct Air Capture (DAC).")
    co2_vent: bool = Field(
        False,
        description="Add option for vent out CO2 from storages to the atmosphere.",
    )
    heat_vent: dict[str, bool] = Field(
        default_factory=lambda: {
            "urban central": True,
            "urban decentral": True,
            "rural": True,
        },
        description="Heat venting by area.",
    )
    marginal_cost_heat_vent: float = Field(
        0.02, description="The marginal cost of heat-venting in all heating systems."
    )

    allam_cycle_gas: bool = Field(
        False,
        description="Add option to include `Allam cycle gas power plants <https://en.wikipedia.org/wiki/Allam_power_cycle>`_.",
    )
    hydrogen_fuel_cell: bool = Field(
        True,
        description="Add option to include hydrogen fuel cell for re-electrification. Assuming OCGT technology costs.",
    )
    hydrogen_turbine: bool = Field(
        True,
        description="Add option to include hydrogen turbine for re-electrification. Assuming OCGT technology costs.",
    )
    SMR: bool = Field(
        True,
        description="Add option for transforming natural gas into hydrogen and CO2 using Steam Methane Reforming (SMR).",
    )
    SMR_cc: bool = Field(
        True,
        description="Add option for transforming natural gas into hydrogen and CO2 using Steam Methane Reforming (SMR) and Carbon Capture (CC).",
    )

    regional_oil_demand: bool = Field(
        True,
        description="Spatially resolve oil demand. Set to true if regional CO2 constraints needed.",
    )
    regional_coal_demand: bool = Field(False, description="Regional coal demand.")

    regional_co2_sequestration_potential: dict[str, Any] = Field(
        default_factory=lambda: {
            "enable": True,
            "attribute": [
                "conservative estimate Mt",
                "conservative estimate GAS Mt",
                "conservative estimate OIL Mt",
                "conservative estimate aquifer Mt",
            ],
            "include_onshore": False,
            "min_size": 3,
            "max_size": 25,
            "years_of_storage": 25,
        },
        description="Add option for regionally-resolved geological carbon dioxide sequestration potentials based on `CO2StoP <https://setis.ec.europa.eu/european-co2-storage-database_en>`_.",
    )
    co2_sequestration_potential: dict[int, float] = Field(
        default_factory=lambda: {
            2020: 0,
            2025: 0,
            2030: 40,
            2035: 100,
            2040: 180,
            2045: 250,
            2050: 250,
        },
        description="The potential of sequestering CO2 in Europe per year and investment period.",
    )
    co2_sequestration_cost: float = Field(
        30, description="The cost of sequestering a ton of CO2 (currency/tCO2)."
    )
    co2_sequestration_lifetime: int = Field(
        50, description="The lifetime of a CO2 sequestration site (years)."
    )
    co2_spatial: bool = Field(
        True,
        description="Add option to spatially resolve carrier representing stored carbon dioxide. This allows for more detailed modelling of CCUTS, e.g. regarding the capturing of industrial process emissions, usage as feedstock for electrofuels, transport of carbon dioxide, and geological sequestration sites.",
    )
    co2_network: bool = Field(
        True,
        description="Add option for planning a new carbon dioxide transmission network.",
    )
    co2_network_cost_factor: float = Field(
        1,
        description="The cost factor for the capital cost of the carbon dioxide transmission network.",
    )
    cc_fraction: float = Field(
        0.9,
        description="The default fraction of CO2 captured with post-combustion capture.",
    )

    hydrogen_underground_storage: bool = Field(
        True,
        description="Add options for storing hydrogen underground. Storage potential depends regionally.",
    )
    hydrogen_underground_storage_locations: list[str] = Field(
        default_factory=lambda: ["onshore", "nearshore"],
        description="The location where hydrogen underground storage can be located. Onshore, nearshore, offshore means it must be located more than 50 km away from the sea, within 50 km of the sea, or within the sea itself respectively.",
    )

    methanol: _MethanolConfig = Field(
        default_factory=_MethanolConfig, description="Methanol configuration."
    )

    ammonia: bool | str = Field(
        True,
        description='Add ammonia as a carrier. It can be either true (copperplated NH3), false (no NH3 carrier) or "regional" (regionalised NH3 without network).',
    )
    min_part_load_electrolysis: float = Field(
        0, description="The minimum unit dispatch (`p_min_pu`) for electrolysis."
    )
    min_part_load_fischer_tropsch: float = Field(
        0.5,
        description="The minimum unit dispatch (`p_min_pu`) for the Fischer-Tropsch process.",
    )
    min_part_load_methanolisation: float = Field(
        0.3,
        description="The minimum unit dispatch (`p_min_pu`) for the methanolisation process.",
    )
    min_part_load_methanation: float = Field(
        0.3, description="Minimum part load methanation."
    )

    use_fischer_tropsch_waste_heat: float = Field(
        0.25,
        description="Add option for using waste heat of Fischer Tropsch in district heating networks.",
    )
    use_haber_bosch_waste_heat: float = Field(
        0.25, description="Use Haber-Bosch waste heat."
    )
    use_methanolisation_waste_heat: float = Field(
        0.25, description="Use methanolisation waste heat."
    )
    use_methanation_waste_heat: float = Field(
        0.25, description="Use methanation waste heat."
    )
    use_fuel_cell_waste_heat: float = Field(
        1,
        description="Add option for using waste heat of fuel cells in district heating networks.",
    )
    use_electrolysis_waste_heat: float = Field(
        0.25,
        description="Add option for using waste heat of electrolysis in district heating networks.",
    )

    electricity_transmission_grid: bool = Field(
        True,
        description="Switch for enabling/disabling the electricity transmission grid.",
    )
    electricity_distribution_grid: bool = Field(
        True,
        description="Add a simplified representation of the exchange capacity between transmission and distribution grid level through a link.",
    )
    electricity_distribution_grid_cost_factor: float = Field(
        1.0,
        description="Multiplies the investment cost of the electricity distribution grid.",
    )
    electricity_grid_connection: bool = Field(
        True,
        description="Add the cost of electricity grid connection for onshore wind and solar.",
    )

    transmission_efficiency: _TransmissionEfficiencyConfig = Field(
        default_factory=_TransmissionEfficiencyConfig,
        description="Transmission efficiency configuration.",
    )

    H2_network: bool = Field(True, description="Add option for new hydrogen pipelines.")
    gas_network: bool = Field(
        True,
        description="Add existing natural gas infrastructure, incl. LNG terminals, production and entry-points. The existing gas network is added with a lossless transport model. A length-weighted `k-edge augmentation algorithm <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation>`_ can be run to add new candidate gas pipelines such that all regions of the model can be connected to the gas network. When activated, all the gas demands are regionally disaggregated as well.",
    )
    H2_retrofit: bool = Field(
        False,
        description="Add option for retrofiting existing pipelines to transport hydrogen.",
    )
    H2_retrofit_capacity_per_CH4: float = Field(
        0.6,
        description="The ratio for H2 capacity per original CH4 capacity of retrofitted pipelines. The `European Hydrogen Backbone (April, 2020) p.15 <https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf>`_ 60% of original natural gas capacity could be used in cost-optimal case as H2 capacity.",
    )
    gas_network_connectivity_upgrade: float = Field(
        1,
        description="The number of desired edge connectivity (k) in the length-weighted `k-edge augmentation algorithm <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation>`_ used for the gas network.",
    )
    gas_distribution_grid: bool = Field(
        True, description="Add a gas distribution grid."
    )
    gas_distribution_grid_cost_factor: float = Field(
        1.0,
        description="Multiplier for the investment cost of the gas distribution grid.",
    )

    biomass_spatial: bool = Field(
        True, description="Add option for resolving biomass demand regionally."
    )
    biomass_transport: bool = Field(
        False, description="Add option for transporting solid biomass between nodes."
    )
    biogas_upgrading: bool = Field(True, description="Biogas upgrading.")
    biogas_upgrading_cc: bool = Field(
        False, description="Add option to capture CO2 from biomass upgrading."
    )

    conventional_generation: dict[str, str] = Field(
        default_factory=lambda: {"OCGT": "gas", "CCGT": "gas"},
        description="Add a more detailed description of conventional carriers. Any power generation requires the consumption of fuel from nodes representing that fuel.",
    )

    biomass_to_liquid: bool = Field(
        True,
        description="Add option for transforming solid biomass into liquid fuel with the same properties as oil.",
    )
    biomass_to_liquid_cc: bool = Field(
        False,
        description="Add option for transforming solid biomass into liquid fuel with the same properties as oil with carbon capture.",
    )
    electrobiofuels: bool = Field(True, description="Electrobiofuels.")
    biosng: bool = Field(
        False,
        description="Add option for transforming solid biomass into synthesis gas with the same properties as natural gas.",
    )
    biosng_cc: bool = Field(
        False,
        description="Add option for transforming solid biomass into synthesis gas with the same properties as natural gas with carbon capture.",
    )
    bioH2: bool = Field(
        False,
        description="Add option for transforming solid biomass into hydrogen with carbon capture.",
    )
    municipal_solid_waste: bool = Field(
        False, description="Add option for municipal solid waste."
    )

    limit_max_growth: _LimitMaxGrowthConfig = Field(
        default_factory=_LimitMaxGrowthConfig,
        description="Limit max growth configuration.",
    )
    enhanced_geothermal: _EnhancedGeothermalConfig = Field(
        default_factory=_EnhancedGeothermalConfig,
        description="Enhanced geothermal configuration.",
    )
    solid_biomass_import: _SolidBiomassImportConfig = Field(
        default_factory=_SolidBiomassImportConfig,
        description="Solid biomass import configuration.",
    )
    imports: _ImportsConfig = Field(
        default_factory=_ImportsConfig, description="Imports configuration."
    )
