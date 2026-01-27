# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Electricity configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#electricity
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from scripts.lib.validation.config._base import ConfigModel


class _OperationalReserveConfig(ConfigModel):
    """Configuration for `electricity.operational_reserve` settings."""

    activate: bool = Field(
        False,
        description="Whether to take operational reserve requirements into account during optimisation.",
    )
    epsilon_load: float = Field(
        0.02,
        description="share of total load.",
    )
    epsilon_vres: float = Field(
        0.02,
        description="share of total renewable supply.",
    )
    contingency: float = Field(
        4000,
        description="Fixed reserve capacity (MW).",
    )


class _MaxHoursConfig(BaseModel):
    """Configuration for `electricity.max_hours` settings."""

    battery: float = Field(
        6,
        description="Maximum state of charge capacity of the battery in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )
    li_ion: float = Field(
        6,
        description="Maximum state of charge capacity of the lithium-ion storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
        alias="li-ion",
    )
    lfp: float = Field(
        6,
        description="Maximum state of charge capacity of the lithium-ion-LFP storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )
    vanadium: float = Field(
        10,
        description="Maximum state of charge capacity of the vanadium-redox-flow storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )
    lair: float = Field(
        12,
        description="Maximum state of charge capacity of the liquid-air storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )
    pair: float = Field(
        24,
        description="Maximum state of charge capacity of the compressed-air-adiabatic storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )
    iron_air: float = Field(
        100,
        description="Maximum state of charge capacity of the iron-air storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
        alias="iron-air",
    )
    H2: float = Field(
        168,
        description="Maximum state of charge capacity of the hydrogen storage in terms of hours at full output capacity `p_nom`. Cf. `PyPSA documentation <https://pypsa.readthedocs.io/en/latest/components.html#storage-unit>`_.",
    )


class _ExtendableCarriersConfig(BaseModel):
    """Configuration for `electricity.extendable_carriers` settings."""

    Generator: list[str] = Field(
        default_factory=lambda: [
            "solar",
            "solar-hsat",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "offwind-float",
            "OCGT",
            "CCGT",
        ],
        description="Defines existing or non-existing conventional and renewable power plants to be extendable during the optimization. Conventional generators can only be built/expanded where already existent today. If a listed conventional carrier is not included in the `conventional_carriers` list, the lower limit of the capacity expansion is set to 0.",
    )
    StorageUnit: list[str] = Field(
        default_factory=list,
        description="Adds extendable storage units at every node/bus after clustering without capacity limits and with zero initial capacity. Supported technologies include battery, H2, li-ion, vanadium, lfp, lair, pair, and iron-air.",
    )
    Store: list[str] = Field(
        default_factory=lambda: ["battery", "H2"],
        description="Adds extendable storage units at every node/bus after clustering without capacity limits and with zero initial capacity. Supported technologies include battery, H2, li-ion, vanadium, lfp, lair, pair, and iron-air.",
    )
    Link: list[str] = Field(
        default_factory=list,
        description="Adds extendable links (H2 pipelines only) at every connection where there are lines or HVDC links without capacity limits and with zero initial capacity. Hydrogen pipelines require hydrogen storage to be modelled as `Store`.",
    )


class _TechnologyMappingConfig(BaseModel):
    """Configuration for `electricity.estimate_renewable_capacities.technology_mapping` settings."""

    Offshore: str = Field(
        "offwind-ac",
        description="PyPSA-Eur carrier that is considered for existing offshore wind technology (IRENA, GEM).",
    )
    Onshore: str = Field(
        "onwind",
        description="PyPSA-Eur carrier that is considered for existing onshore wind capacities (IRENA, GEM).",
    )
    PV: str = Field(
        "solar",
        description="PyPSA-Eur carrier that is considered for existing solar PV capacities (IRENA, GEM).",
    )


class _EstimateRenewableCapacitiesConfig(BaseModel):
    """Configuration for `electricity.estimate_renewable_capacities` settings."""

    enable: bool = Field(
        True,
        description="Activate routine to estimate renewable capacities in rule `add_electricity`. This option should not be used in combination with pathway planning `foresight: myopic` or `foresight: perfect` as renewable capacities are added differently in `add_existing_baseyear`.",
    )
    from_gem: bool = Field(
        True,
        description="Add renewable capacities from `Global Energy Monitor's Global Solar Power Tracker <https://globalenergymonitor.org/projects/global-solar-power-tracker/>`_ and `Global Energy Monitor's Global Wind Power Tracker <https://globalenergymonitor.org/projects/global-wind-power-tracker/>`_.",
    )
    year: int = Field(
        2020,
        description="Renewable capacities are based on existing capacities reported by IRENA (IRENASTAT) for the specified year.",
    )
    expansion_limit: float | bool = Field(
        False,
        description="Artificially limit maximum IRENA capacities to a factor. For example, an `expansion_limit: 1.1` means 110% of capacities. If false are chosen, the estimated renewable potentials determine by the workflow are used.",
    )
    technology_mapping: _TechnologyMappingConfig = Field(
        default_factory=_TechnologyMappingConfig,
        description="Mapping between PyPSA-Eur and powerplantmatching technology names.",
    )


class _AutarkyConfig(BaseModel):
    """Configuration for `electricity.autarky` settings."""

    enable: bool = Field(
        False,
        description="Require each node to be autarkic by removing all lines and links.",
    )
    by_country: bool = Field(
        False,
        description="Require each country to be autarkic by removing all cross-border lines and links. `electricity: autarky` must be enabled.",
    )


class ElectricityConfig(BaseModel):
    """Configuration for `electricity` settings."""

    voltages: list[float] = Field(
        default_factory=lambda: [220.0, 300.0, 330.0, 380.0, 400.0, 500.0, 750.0],
        description="Voltage levels to consider.",
    )
    base_network: Literal["entsoegridkit", "osm", "tyndp"] = Field(
        "osm",
        description="Specify the underlying base network, i.e. GridKit (based on ENTSO-E web map extract), OpenStreetMap (OSM), or TYNDP.",
    )
    gaslimit_enable: bool = Field(
        False,
        description="Add an overall absolute gas limit configured in `electricity: gaslimit`.",
    )
    gaslimit: float | bool = Field(
        False,
        description="Global gas usage limit.",
    )
    co2limit_enable: bool = Field(
        False,
        description="Add an overall absolute carbon-dioxide emissions limit configured in `electricity: co2limit` in `prepare_network`. **Warning:** This option should currently only be used with electricity-only networks, not for sector-coupled networks.",
    )
    co2limit: float = Field(
        7.75e7,
        description="Cap on total annual system carbon dioxide emissions.",
    )
    co2base: float = Field(
        1.487e9,
        description="Reference value of total annual system carbon dioxide emissions if relative emission reduction target is specified in `{opts}` wildcard.",
    )
    operational_reserve: _OperationalReserveConfig = Field(
        default_factory=_OperationalReserveConfig,
        description="Settings for reserve requirements following `GenX <https://genxproject.github.io/GenX/dev/core/#Reserves>`_.",
    )
    max_hours: _MaxHoursConfig = Field(
        default_factory=_MaxHoursConfig,
        description="Maximum storage hours configuration.",
    )
    extendable_carriers: _ExtendableCarriersConfig = Field(
        default_factory=_ExtendableCarriersConfig,
        description="Defines which carriers are extendable during optimization.",
    )
    powerplants_filter: str | bool = Field(
        "(DateOut >= 2024 or DateOut != DateOut) and not (Country == 'Germany' and Fueltype == 'Nuclear')",
        description="Filter query for the default powerplant database.",
    )
    custom_powerplants: str | bool = Field(
        False,
        description="Filter query for the custom powerplant database.",
    )
    everywhere_powerplants: list[str] = Field(
        default_factory=list,
        description="List of conventional power plants to add to every node in the model with zero initial capacity. To be used in combination with `extendable_carriers` to allow for building conventional powerplants irrespective of existing locations.",
    )
    conventional_carriers: list[str] = Field(
        default_factory=lambda: [
            "nuclear",
            "oil",
            "OCGT",
            "CCGT",
            "coal",
            "lignite",
            "geothermal",
            "biomass",
        ],
        description="List of conventional power plants to include in the model from `resources/powerplants_s_{clusters}.csv`. If an included carrier is also listed in `extendable_carriers`, the capacity is taken as a lower bound.",
    )
    renewable_carriers: list[str] = Field(
        default_factory=lambda: [
            "solar",
            "solar-hsat",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "offwind-float",
            "hydro",
        ],
        description="List of renewable generators to include in the model.",
    )
    estimate_renewable_capacities: _EstimateRenewableCapacitiesConfig = Field(
        default_factory=_EstimateRenewableCapacitiesConfig,
        description="Configuration for estimating renewable capacities.",
    )
    autarky: _AutarkyConfig = Field(
        default_factory=_AutarkyConfig,
        description="Autarky configuration.",
    )
    transmission_limit: str = Field(
        "vopt",
        description="Limit on transmission expansion. The first part can be `v` (for setting a limit on line volume) or `c` (for setting a limit on line cost). The second part can be `opt` or a float bigger than one (e.g. 1.25). If `opt` is chosen line expansion is optimised according to its capital cost (where the choice `v` only considers overhead costs for HVDC transmission lines, while `c` uses more accurate costs distinguishing between overhead and underwater sections and including inverter pairs). The setting `v1.25` will limit the total volume of line expansion to 25% of currently installed capacities weighted by individual line lengths. The setting `c1.25` will allow to build a transmission network that costs no more than 25 % more than the current system.",
    )

    model_config = ConfigDict(populate_by_name=True)
