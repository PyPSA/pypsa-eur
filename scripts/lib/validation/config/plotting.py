# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Plotting configuration block.

See # docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#plotting
"""

from typing import Literal

from pydantic import ConfigDict, Field

from scripts.lib.validation.config._base import ConfigModel


class _MapConfig(ConfigModel):
    """Configuration for `plotting.map` settings."""

    boundaries: list[float] | None = Field(
        [-11, 30, 34, 71],
        description="Bounding box of static maps as ``[x_min, x_max, y_min, y_max]`` in longitude/latitude degrees. If ``null``, the boundaries are derived automatically from the plotted regions.",
    )
    geomap_colors: dict[str, str] = Field(
        default_factory=lambda: {"ocean": "white", "land": "white"},
        description="Colors used for cartopy geographic map features, e.g. ``ocean``, ``land``, ``border`` or ``coastline``.",
    )


class _ProjectionConfig(ConfigModel):
    """Configuration for `plotting.projection` settings."""

    model_config = ConfigDict(extra="allow")

    name: str = Field(
        "EqualEarth",
        description="Name of the `cartopy <https://scitools.org.uk/cartopy/docs/latest/reference/projections.html>`_ CRS class used for the map projection, e.g. ``EqualEarth`` or ``LambertConformal``. Additional keys are passed on as keyword arguments to the projection class, e.g. ``central_longitude``, ``central_latitude`` or ``standard_parallels`` for ``LambertConformal``.",
    )


class _EuNodeLocationConfig(ConfigModel):
    """Configuration for `plotting.eu_node_location` settings."""

    x: float = Field(
        -5.5,
        description="Longitude of the auxiliary EU-wide node used for continent-level summaries.",
    )
    y: float = Field(
        46.0,
        description="Latitude of the auxiliary EU-wide node used for continent-level summaries.",
    )


class _BalanceTimeseriesConfig(ConfigModel):
    """Configuration for `plotting.balance_timeseries` settings."""

    max_threshold: float = Field(
        5,
        description="Maximum absolute value (in GW) used to determine the y-axis limits of balance timeseries plots.",
    )
    mean_threshold: float = Field(
        1,
        description="Mean absolute value (in GW) below which a carrier is excluded from balance timeseries plots.",
    )
    monthly: bool = Field(
        True, description="Enable generation of monthly balance timeseries plots."
    )
    monthly_resolution: str | None = Field(
        None,
        description="Resampling resolution (pandas offset alias) for monthly balance timeseries plots. ``null`` uses the native snapshot resolution.",
    )
    annual: bool = Field(
        True, description="Enable generation of annual balance timeseries plots."
    )
    annual_resolution: str = Field(
        "D",
        description="Resampling resolution (pandas offset alias, e.g. ``D`` for daily) for annual balance timeseries plots.",
    )
    carriers: list[str] = Field(
        [
            "H2",
            "NH3",
            "gas",
            "methanol",
            "oil",
            "solid biomass",
            "biogas",
            "co2 stored",
            "co2",
        ],
        description="Bus carriers for which balance timeseries plots are generated.",
    )
    carrier_groups: dict[str, list[str]] = Field(
        default_factory=lambda: {
            "electricity": ["AC", "low voltage"],
            "heat": [
                "urban central heat",
                "urban decentral heat",
                "rural heat",
                "residential urban decentral heat",
                "residential rural heat",
                "services urban decentral heat",
                "services rural heat",
            ],
        },
        description="Named groups of bus carriers that are aggregated into a single balance timeseries plot.",
    )


class _InteractiveBusBalanceConfig(ConfigModel):
    """Configuration for `plotting.interactive_bus_balance` settings."""

    bus_name_pattern: str | None = Field(
        "None",
        description="Shell-style glob pattern (see Python's ``fnmatch``) used to select buses for the interactive bus balance plot, e.g. ``DE*`` for German buses. The default value matches no bus name literally, so no plots are generated unless overridden. Set to ``NONE_BY_DEFAULT`` to explicitly disable, or ``null``/empty to include all buses.",
    )


class _HeatmapTimeseriesConfig(ConfigModel):
    """Configuration for `plotting.heatmap_timeseries` settings."""

    marginal_price: list[str] = Field(
        [
            "AC",
            "H2",
            "NH3",
            "gas",
            "methanol",
            "oil",
            "co2 stored",
            "urban central heat",
        ],
        description="Bus carriers for which marginal price heatmaps are generated.",
    )
    utilisation_rate: list[str] = Field(
        [
            "solar",
            "solar rooftop",
            "solar-hsat",
            "onwind",
            "offwind-dc",
            "offwind-ac",
            "offwind-float",
            "ror",
            "hydro",
            "PHS",
            "battery charger",
            "battery discharger",
            "H2 Electrolysis",
            "Fischer-Tropsch",
            "methanolisation",
            "Sabatier",
            "OCGT",
            "H2 Fuel Cell",
            "urban central CHP",
            "urban central CHP CC",
            "urban central solid biomass CHP",
            "urban central solid biomass CHP CC",
            "rural gas boiler",
            "urban central air heat pump",
            "DAC",
        ],
        description="Generator/storage carriers for which utilisation-rate heatmaps are generated.",
    )
    soc: list[str] = Field(
        [
            "battery",
            "H2 Store",
            "co2 stored",
            "gas",
            "methanol",
            "oil",
            "urban central water tanks",
        ],
        description="Storage carriers for which state-of-charge heatmaps are generated.",
    )


class _BalanceMapCarrierConfig(ConfigModel):
    """Style settings for a single carrier of `plotting.balance_map`."""

    cmap: str = Field(
        description="Matplotlib colormap used to color regions by nodal price/marginal value."
    )
    vmin: float | None = Field(
        None,
        description="Lower bound of the region color scale. ``null`` scales automatically to the data.",
    )
    vmax: float | None = Field(
        None,
        description="Upper bound of the region color scale. ``null`` scales automatically to the data.",
    )
    region_unit: str = Field(
        description="Unit label shown in the region color bar legend."
    )
    branch_color: str = Field(
        description="Color used for transmission branches (lines/links) carrying this carrier."
    )
    unit: str = Field(description="Unit label shown in the bus/branch size legend.")
    unit_conversion: float = Field(
        description="Factor used to convert the carrier's native unit (MWh) into ``unit``."
    )
    bus_factor: float = Field(
        description="Scaling factor applied to bus pie/marker sizes."
    )
    branch_factor: float = Field(
        description="Scaling factor applied to branch line widths."
    )
    flow_factor: float = Field(
        description="Scaling factor applied to flow arrow sizes."
    )
    bus_sizes: list[float] = Field(
        description="Reference values (in ``unit``) shown in the bus size legend."
    )
    branch_sizes: list[float] | None = Field(
        None,
        description="Reference values (in ``unit``) shown in the branch size legend. ``null`` omits the branch size legend.",
    )


class _BalanceMapConfig(ConfigModel):
    """
    Configuration for `plotting.balance_map` settings.

    Besides ``bus_carriers``, arbitrary additional keys are allowed, one per bus carrier,
    each configuring the static balance map style for that carrier (see
    :class:`_BalanceMapCarrierConfig`). They are looked up dynamically by carrier name.
    """

    model_config = ConfigDict(extra="allow")
    __pydantic_extra__: dict[str, _BalanceMapCarrierConfig]

    bus_carriers: list[str] = Field(
        ["AC", "co2_stored", "gas", "H2", "methanol", "oil", "urban_central_heat"],
        description="Bus carriers for which a static balance map is generated (see `plot_balance_map`).",
    )


class _BalanceMapInteractiveCarrierConfig(ConfigModel):
    """Style settings for a single carrier of `plotting.balance_map_interactive`."""

    cmap: str = Field(
        description="Matplotlib colormap used to color regions by nodal price/marginal value."
    )
    region_unit: str = Field(
        description="Unit label shown in the region color bar legend."
    )
    vmin: float | None = Field(
        None,
        description="Lower bound of the region color scale. ``null`` scales automatically to the data.",
    )
    vmax: float | None = Field(
        None,
        description="Upper bound of the region color scale. ``null`` scales automatically to the data.",
    )
    region_alpha: float = Field(description="Opacity of the region choropleth layer.")
    unit_conversion: float = Field(
        description="Factor used to convert the carrier's native unit (MWh) into the legend unit."
    )
    branch_color: str = Field(
        description="Color used for transmission branches (lines/links) carrying this carrier."
    )
    branch_width_max: float = Field(description="Maximum branch line width in pixels.")
    bus_size_max: float = Field(description="Maximum bus marker size in pixels.")
    arrow_size_factor: float = Field(
        description="Scaling factor applied to flow arrow sizes."
    )
    map_style: str = Field(
        description="Base map tile style passed to the plotting backend, e.g. ``road`` or ``satellite``."
    )
    tooltip: bool = Field(
        description="Enable hover tooltips showing bus/branch values."
    )


class _BalanceMapInteractiveConfig(ConfigModel):
    """
    Configuration for `plotting.balance_map_interactive` settings.

    Besides ``bus_carriers``, arbitrary additional keys are allowed, one per bus carrier,
    each configuring the interactive balance map style for that carrier (see
    :class:`_BalanceMapInteractiveCarrierConfig`). They are looked up dynamically by carrier
    name.
    """

    model_config = ConfigDict(extra="allow")
    __pydantic_extra__: dict[str, _BalanceMapInteractiveCarrierConfig]

    bus_carriers: list[str] = Field(
        ["AC", "co2_stored", "gas", "H2", "methanol", "oil", "urban_central_heat"],
        description="Bus carriers for which an interactive balance map is generated (see `plot_balance_map_interactive`).",
    )


class _HeatSourceMapConfig(ConfigModel):
    """Configuration for `plotting.heat_source_map` settings."""

    temperature_cmap: str = Field(
        "Reds", description="Matplotlib colormap for the heat source temperature map."
    )
    energy_cmap: str = Field(
        "Oranges",
        description="Matplotlib colormap for the heat source energy potential map.",
    )


def _default_balance_map() -> _BalanceMapConfig:
    """Default value for `plotting.balance_map`."""
    return _BalanceMapConfig(
        bus_carriers=[
            "AC",
            "co2_stored",
            "gas",
            "H2",
            "methanol",
            "oil",
            "urban_central_heat",
        ],
        AC=_BalanceMapCarrierConfig(
            cmap="Greens",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="darkseagreen",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.002,
            branch_factor=0.01,
            flow_factor=100,
            bus_sizes=[200, 100],
            branch_sizes=[100, 20],
        ),
        biogas=_BalanceMapCarrierConfig(
            cmap="Greens",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="darkseagreen",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.1,
            branch_factor=0.1,
            flow_factor=100,
            bus_sizes=[100, 50],
            branch_sizes=None,
        ),
        co2_stored=_BalanceMapCarrierConfig(
            cmap="Purples",
            vmin=None,
            vmax=None,
            region_unit="€/t_${CO_2}$",
            branch_color="orange",
            unit="Mt",
            unit_conversion=1000000,
            bus_factor=0.015,
            branch_factor=0.4,
            flow_factor=120,
            bus_sizes=[50, 10],
            branch_sizes=[5, 2],
        ),
        gas=_BalanceMapCarrierConfig(
            cmap="Oranges",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="darkred",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.002,
            branch_factor=0.05,
            flow_factor=60,
            bus_sizes=[200, 100],
            branch_sizes=[100, 50],
        ),
        H2=_BalanceMapCarrierConfig(
            cmap="Blues",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="pink",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.002,
            branch_factor=0.03,
            flow_factor=8,
            bus_sizes=[50, 25],
            branch_sizes=[40, 20],
        ),
        methanol=_BalanceMapCarrierConfig(
            cmap="Greens",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="yellow",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.005,
            branch_factor=0.1,
            flow_factor=100,
            bus_sizes=[20, 10],
            branch_sizes=None,
        ),
        oil=_BalanceMapCarrierConfig(
            cmap="Greys",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            branch_color="black",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.002,
            branch_factor=0.01,
            flow_factor=100,
            bus_sizes=[200, 100],
            branch_sizes=None,
        ),
        solid_biomass=_BalanceMapCarrierConfig(
            cmap="Greens",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="darkseagreen",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.01,
            branch_factor=0.1,
            flow_factor=100,
            bus_sizes=[100, 50],
            branch_sizes=None,
        ),
        urban_central_heat=_BalanceMapCarrierConfig(
            cmap="Oranges",
            vmin=None,
            vmax=None,
            region_unit="€/MWh",
            branch_color="darkred",
            unit="TWh",
            unit_conversion=1000000,
            bus_factor=0.005,
            branch_factor=0.1,
            flow_factor=100,
            bus_sizes=[300, 100],
            branch_sizes=None,
        ),
    )


def _default_balance_map_interactive() -> _BalanceMapInteractiveConfig:
    """Default value for `plotting.balance_map_interactive`."""
    return _BalanceMapInteractiveConfig(
        bus_carriers=[
            "AC",
            "co2_stored",
            "gas",
            "H2",
            "methanol",
            "oil",
            "urban_central_heat",
        ],
        AC=_BalanceMapInteractiveCarrierConfig(
            cmap="Greens",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="darkseagreen",
            branch_width_max=20,
            bus_size_max=15000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        co2_stored=_BalanceMapInteractiveCarrierConfig(
            cmap="Purples",
            region_unit="€/t CO2",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="orange",
            branch_width_max=20,
            bus_size_max=15000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        gas=_BalanceMapInteractiveCarrierConfig(
            cmap="Oranges",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="darkred",
            branch_width_max=20,
            bus_size_max=15000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        H2=_BalanceMapInteractiveCarrierConfig(
            cmap="Blues",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="pink",
            branch_width_max=45,
            bus_size_max=20000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        methanol=_BalanceMapInteractiveCarrierConfig(
            cmap="Greens",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="yellow",
            branch_width_max=45,
            bus_size_max=20000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        oil=_BalanceMapInteractiveCarrierConfig(
            cmap="Greys",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="black",
            branch_width_max=45,
            bus_size_max=20000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        solid_biomass=_BalanceMapInteractiveCarrierConfig(
            cmap="Greens",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="darkseagreen",
            branch_width_max=45,
            bus_size_max=20000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
        urban_central_heat=_BalanceMapInteractiveCarrierConfig(
            cmap="Oranges",
            region_unit="€/MWh",
            vmin=None,
            vmax=None,
            region_alpha=0.8,
            unit_conversion=1000000,
            branch_color="darkred",
            branch_width_max=45,
            bus_size_max=20000,
            arrow_size_factor=2,
            map_style="road",
            tooltip=True,
        ),
    )


_NICE_NAMES: dict[str, str] = {
    "OCGT": "Open-Cycle Gas",
    "CCGT": "Combined-Cycle Gas",
    "offwind-ac": "Offshore Wind (AC)",
    "offwind-dc": "Offshore Wind (DC)",
    "offwind-float": "Offshore Wind (Floating)",
    "onwind": "Onshore Wind",
    "solar": "Solar",
    "PHS": "Pumped Hydro Storage",
    "hydro": "Reservoir & Dam",
    "battery": "Battery Storage",
    "H2": "Hydrogen Storage",
    "li-ion": "Lithium-Ion Storage",
    "lfp": "Lithium-Ion-LFP Storage",
    "vanadium": "Vanadium-Redox-Flow Storage",
    "lair": "Liquid-Air Storage",
    "pair": "Compressed-Air-Adiabatic Storage",
    "iron-air": "Iron-Air Storage",
    "lines": "Transmission Lines",
    "ror": "Run of River",
    "load": "Load Shedding",
    "ac": "AC",
    "dc": "DC",
}


_TECH_COLORS: dict[str, str] = {
    # wind
    "onwind": "#235ebc",
    "onshore wind": "#235ebc",
    "offwind": "#6895dd",
    "offshore wind": "#6895dd",
    "offwind-ac": "#6895dd",
    "offshore wind (AC)": "#6895dd",
    "offshore wind ac": "#6895dd",
    "offwind-dc": "#74c6f2",
    "offshore wind (DC)": "#74c6f2",
    "offshore wind dc": "#74c6f2",
    "offwind-float": "#b5e2fa",
    "offshore wind (Float)": "#b5e2fa",
    "offshore wind float": "#b5e2fa",
    # water
    "hydro": "#298c81",
    "hydro reservoir": "#298c81",
    "ror": "#3dbfb0",
    "run of river": "#3dbfb0",
    "hydroelectricity": "#298c81",
    "PHS": "#51dbcc",
    "hydro+PHS": "#08ad97",
    # solar
    "solar": "#f9d002",
    "solar PV": "#f9d002",
    "solar-hsat": "#fdb915",
    "solar thermal": "#ffbf2b",
    "residential rural solar thermal": "#f1c069",
    "services rural solar thermal": "#eabf61",
    "residential urban decentral solar thermal": "#e5bc5a",
    "services urban decentral solar thermal": "#dfb953",
    "urban central solar thermal": "#d7b24c",
    "solar rooftop": "#ffea80",
    # gas
    "OCGT": "#e0986c",
    "OCGT marginal": "#e0986c",
    "OCGT-heat": "#e0986c",
    "gas boiler": "#db6a25",
    "gas boilers": "#db6a25",
    "gas boiler marginal": "#db6a25",
    "residential rural gas boiler": "#d4722e",
    "residential urban decentral gas boiler": "#cb7a36",
    "services rural gas boiler": "#c4813f",
    "services urban decentral gas boiler": "#ba8947",
    "urban central gas boiler": "#b0904f",
    "gas": "#e05b09",
    "fossil gas": "#e05b09",
    "natural gas": "#e05b09",
    "biogas to gas": "#e36311",
    "biogas to gas CC": "#e51245",
    "CCGT": "#a85522",
    "CCGT marginal": "#a85522",
    "allam": "#B98F76",
    "gas for industry co2 to atmosphere": "#692e0a",
    "gas for industry co2 to stored": "#8a3400",
    "gas for industry": "#853403",
    "gas for industry CC": "#692e0a",
    "gas pipeline": "#ebbca0",
    "gas pipeline new": "#a87c62",
    # oil
    "oil": "#c9c9c9",
    "oil primary": "#d2d2d2",
    "oil refining": "#e6e6e6",
    "imported oil": "#a3a3a3",
    "oil boiler": "#adadad",
    "residential rural oil boiler": "#a9a9a9",
    "services rural oil boiler": "#a5a5a5",
    "residential urban decentral oil boiler": "#a1a1a1",
    "urban central oil boiler": "#9d9d9d",
    "services urban decentral oil boiler": "#999999",
    "agriculture machinery oil": "#949494",
    "agriculture machinery electric": "#444578",
    "shipping oil": "#808080",
    "land transport oil": "#afafaf",
    # nuclear
    "Nuclear": "#ff8c00",
    "Nuclear marginal": "#ff8c00",
    "nuclear": "#ff8c00",
    "uranium": "#ff8c00",
    # coal
    "Coal": "#545454",
    "coal": "#545454",
    "Coal marginal": "#545454",
    "coal for industry": "#343434",
    "solid": "#545454",
    "Lignite": "#826837",
    "lignite": "#826837",
    "Lignite marginal": "#826837",
    # biomass
    "biogas": "#e3d37d",
    "biomass": "#baa741",
    "solid biomass": "#baa741",
    "municipal solid waste": "#91ba41",
    "solid biomass import": "#d5ca8d",
    "solid biomass transport": "#baa741",
    "solid biomass for industry": "#7a6d26",
    "solid biomass for industry CC": "#47411c",
    "solid biomass for industry co2 from atmosphere": "#736412",
    "solid biomass for industry co2 to stored": "#47411c",
    "urban central solid biomass CHP": "#9d9042",
    "urban central solid biomass CHP CC": "#6c5d28",
    "biomass boiler": "#8A9A5B",
    "residential rural biomass boiler": "#a1a066",
    "residential urban decentral biomass boiler": "#b0b87b",
    "services rural biomass boiler": "#c6cf98",
    "services urban decentral biomass boiler": "#dde5b5",
    "biomass to liquid": "#32CD32",
    "unsustainable solid biomass": "#998622",
    "unsustainable bioliquids": "#32CD32",
    "electrobiofuels": "red",
    "BioSNG": "#123456",
    "BioSNG CC": "#45233b",
    "solid biomass to hydrogen": "#654321",
    # power transmission
    "lines": "#6c9459",
    "transmission lines": "#6c9459",
    "electricity distribution grid": "#97ad8c",
    "low voltage": "#97ad8c",
    # electricity demand
    "Electric load": "#110d63",
    "electric demand": "#110d63",
    "electricity": "#110d63",
    "industry electricity": "#2d2a66",
    "industry new electricity": "#2d2a66",
    "agriculture electricity": "#494778",
    # battery + EVs
    "battery": "#ace37f",
    "battery storage": "#ace37f",
    "battery charger": "#88a75b",
    "battery discharger": "#5d4e29",
    "home battery": "#80c944",
    "home battery storage": "#80c944",
    "home battery charger": "#5e8032",
    "home battery discharger": "#3c5221",
    "BEV charger": "#baf238",
    "V2G": "#e5ffa8",
    "land transport EV": "#baf238",
    "land transport demand": "#38baf2",
    "EV battery": "#baf238",
    # all battery variation:
    "li-ion": "#ace37f",
    "li-ion charger": "#ace37f",
    "li-ion discharger": "#ace37f",
    "lfp": "#ace37f",
    "lfp charger": "#ace37f",
    "lfp discharger": "#ace37f",
    "vanadium": "#9B111E",
    "vanadium charger": "#9B111E",
    "vanadium discharger": "#9B111E",
    "lair": "#87CEEB",
    "lair charger": "#87CEEB",
    "lair discharger": "#87CEEB",
    "pair": "#003366",
    "pair charger": "#003366",
    "pair discharger": "#003366",
    "iron-air": "#edba1c",
    "iron-air charger": "#edba1c",
    "iron-air discharger": "#edba1c",
    # hot water storage
    "water tanks": "#e69487",
    "residential rural water tanks": "#f7b7a3",
    "services rural water tanks": "#f3afa3",
    "residential urban decentral water tanks": "#f2b2a3",
    "services urban decentral water tanks": "#f1b4a4",
    "urban central water tanks": "#e9977d",
    "hot water storage": "#e69487",
    "hot water charging": "#e8998b",
    "urban central water tanks charger": "#b57a67",
    "residential rural water tanks charger": "#b4887c",
    "residential urban decentral water tanks charger": "#b39995",
    "services rural water tanks charger": "#b3abb0",
    "services urban decentral water tanks charger": "#b3becc",
    "hot water discharging": "#e99c8e",
    "urban central water tanks discharger": "#b9816e",
    "residential rural water tanks discharger": "#ba9685",
    "residential urban decentral water tanks discharger": "#baac9e",
    "services rural water tanks discharger": "#bbc2b8",
    "services urban decentral water tanks discharger": "#bdd8d3",
    "water pits": "#cc826a",
    "water pits charger": "#b36a5e",
    "water pits discharger": "#b37468",
    "urban central water pits": "#d96f4c",
    "urban central water pits charger": "#a85d47",
    "urban central water pits discharger": "#b36452",
    "aquifer thermal energy storage": "#6d00fc",
    "aquifer thermal energy storage charger": "#6d00fc",
    "aquifer thermal energy storage discharger": "#6d00fc",
    # heat demand
    "Heat load": "#cc1f1f",
    "heat": "#cc1f1f",
    "heat vent": "#aa3344",
    "heat demand": "#cc1f1f",
    "rural heat": "#ff5c5c",
    "rural heat dsm": "#ff5c5c",
    "heat dsm": "#ff5c5c",
    "residential rural heat": "#ff7c7c",
    "services rural heat": "#ff9c9c",
    "central heat": "#cc1f1f",
    "urban central heat": "#d15959",
    "urban central heat dsm": "#d15959",
    "urban central heat vent": "#a74747",
    "decentral heat": "#750606",
    "residential urban decentral heat": "#a33c3c",
    "residential urban decentral heat dsm": "#a33c3c",
    "services urban decentral heat": "#cc1f1f",
    "low-temperature heat for industry": "#8f2727",
    "process heat": "#ff0000",
    "agriculture heat": "#d9a5a5",
    # heat supply
    "heat pumps": "#2fb537",
    "heat pump": "#2fb537",
    "air heat pump": "#36eb41",
    "residential urban decentral air heat pump": "#48f74f",
    "services urban decentral air heat pump": "#5af95d",
    "services rural air heat pump": "#5af95d",
    "urban central air heat pump": "#6cfb6b",
    "ptes heat pump": "#5dade2",
    "urban central ptes heat pump": "#3498db",
    "urban central geothermal heat pump": "#4f2144",
    "geothermal heat pump": "#4f2144",
    "geothermal heat direct utilisation": "#ba91b1",
    "river_water heat": "#4bb9f2",
    "river_water heat pump": "#4bb9f2",
    "sea_water heat": "#0b222e",
    "sea_water heat pump": "#0b222e",
    "ground heat pump": "#2fb537",
    "residential rural ground heat pump": "#4f2144",
    "residential rural air heat pump": "#48f74f",
    "services rural ground heat pump": "#5af95d",
    "Ambient": "#98eb9d",
    "CHP": "#8a5751",
    "urban central gas CHP": "#8d5e56",
    "CHP CC": "#634643",
    "urban central gas CHP CC": "#6e4e4c",
    "CHP heat": "#8a5751",
    "CHP electric": "#8a5751",
    "district heating": "#e8beac",
    "resistive heater": "#d8f9b8",
    "residential rural resistive heater": "#bef5b5",
    "residential urban decentral resistive heater": "#b2f1a9",
    "services rural resistive heater": "#a5ed9d",
    "services urban decentral resistive heater": "#98e991",
    "urban central resistive heater": "#8cdf85",
    "retrofitting": "#8487e8",
    "building retrofitting": "#8487e8",
    # hydrogen
    "H2 for industry": "#f073da",
    "H2 for shipping": "#ebaee0",
    "H2": "#bf13a0",
    "hydrogen": "#bf13a0",
    "retrofitted H2 boiler": "#e5a0d9",
    "SMR": "#870c71",
    "SMR CC": "#4f1745",
    "H2 liquefaction": "#d647bd",
    "hydrogen storage": "#bf13a0",
    "H2 Store": "#bf13a0",
    "H2 storage": "#bf13a0",
    "land transport fuel cell": "#6b3161",
    "H2 pipeline": "#f081dc",
    "H2 pipeline retrofitted": "#ba99b5",
    "H2 Fuel Cell": "#c251ae",
    "H2 fuel cell": "#c251ae",
    "H2 turbine": "#991f83",
    "H2 Electrolysis": "#ff29d9",
    "H2 electrolysis": "#ff29d9",
    # ammonia
    "NH3": "#46caf0",
    "ammonia": "#46caf0",
    "ammonia store": "#00ace0",
    "ammonia cracker": "#87d0e6",
    "Haber-Bosch": "#076987",
    # syngas
    "Sabatier": "#9850ad",
    "methanation": "#c44ce6",
    "methane": "#c44ce6",
    # synfuels
    "Fischer-Tropsch": "#25c49a",
    "liquid": "#25c49a",
    "kerosene for aviation": "#a1ffe6",
    "naphtha for industry": "#57ebc4",
    "methanol-to-kerosene": "#C98468",
    "methanol-to-olefins/aromatics": "#FFA07A",
    "Methanol steam reforming": "#FFBF00",
    "Methanol steam reforming CC": "#A2EA8A",
    "methanolisation": "#00FFBF",
    "biomass-to-methanol": "#EAD28A",
    "biomass-to-methanol CC": "#EADBAD",
    "allam methanol": "#B98F76",
    "CCGT methanol": "#B98F76",
    "CCGT methanol CC": "#B98F76",
    "OCGT methanol": "#B98F76",
    "methanol": "#FF7B00",
    "methanol transport": "#FF7B00",
    "shipping methanol": "#468c8b",
    "industry methanol": "#468c8b",
    # co2
    "CC": "#f29dae",
    "CCS": "#f29dae",
    "CO2 sequestration": "#f29dae",
    "DAC": "#ff5270",
    "co2 stored": "#f2385a",
    "co2 sequestered": "#f2682f",
    "co2 dense": "#65334d",
    "co2 expansion": "#c6ebbe",
    "co2 compression": "#a9dbb8",
    "co2": "#f29dae",
    "co2 vent": "#ffd4dc",
    "CO2 pipeline": "#f5627f",
    # emissions
    "process emissions CC": "#000000",
    "process emissions": "#222222",
    "process emissions to stored": "#444444",
    "process emissions to atmosphere": "#888888",
    "oil emissions": "#aaaaaa",
    "shipping oil emissions": "#555555",
    "shipping methanol emissions": "#666666",
    "land transport oil emissions": "#777777",
    "agriculture machinery oil emissions": "#333333",
    # other
    "shipping": "#03a2ff",
    "power-to-heat": "#2fb537",
    "power-to-gas": "#c44ce6",
    "power-to-H2": "#ff29d9",
    "power-to-liquid": "#25c49a",
    "gas-to-power/heat": "#ee8340",
    "waste": "#e3d37d",
    "other": "#000000",
    "geothermal": "#ba91b1",
    "geothermal heat": "#ba91b1",
    "geothermal district heat": "#d19D00",
    "geothermal organic rankine cycle": "#ffbf00",
    "AC": "#70af1d",
    "AC-AC": "#70af1d",
    "AC line": "#70af1d",
    "links": "#8a1caf",
    "HVDC links": "#8a1caf",
    "DC": "#8a1caf",
    "DC-DC": "#8a1caf",
    "DC link": "#8a1caf",
    "load": "#dd2e23",
    "waste CHP": "#e3d37d",
    "waste CHP CC": "#e3d3ff",
    "non-sequestered HVC": "#8f79b5",
    "HVC to air": "k",
    "import H2": "#db8ccd",
    "import gas": "#f7a572",
    "import NH3": "#e2ed74",
    "import oil": "#93eda2",
    "import methanol": "#87d0e6",
}


class PlottingConfig(ConfigModel):
    """Configuration for top level `plotting` settings."""

    enable_heat_source_maps: bool = Field(
        False,
        description="Enable generation of heat source temperature and energy potential maps (see `plot_heat_source_map`).",
    )
    heat_source_map: _HeatSourceMapConfig = Field(
        default_factory=_HeatSourceMapConfig,
        description="Colormap settings for heat source temperature/energy maps, used when ``enable_heat_source_maps`` is set.",
    )
    map: _MapConfig = Field(
        default_factory=_MapConfig,
        description="Settings for the geographic extent and base map colors of static maps.",
    )
    projection: _ProjectionConfig = Field(
        default_factory=_ProjectionConfig,
        description="Cartopy map projection settings for static maps.",
    )
    eu_node_location: _EuNodeLocationConfig = Field(
        default_factory=_EuNodeLocationConfig,
        description="Location of the auxiliary EU-wide node used for continent-level summaries.",
    )
    costs_max: float | Literal["auto"] | None = Field(
        1000.0,
        description="Upper limit (in billion EUR/a) for the y-axis of cost bar plots. Set to ``auto`` to scale automatically to the data, or ``null`` for the plotting library default.",
    )
    costs_threshold: float = Field(
        1,
        description="Minimum cost (in billion EUR/a) for a cost component to be included in cost plots.",
    )
    energy_max: float | Literal["auto"] | None = Field(
        20000.0,
        description="Upper limit (in TWh) for the y-axis of energy balance bar plots. Set to ``auto`` to scale automatically to the data, or ``null`` for the plotting library default.",
    )
    energy_min: float | Literal["auto"] | None = Field(
        -20000.0,
        description="Lower limit (in TWh) for the y-axis of energy balance bar plots. Set to ``auto`` to scale automatically to the data, or ``null`` for the plotting library default.",
    )
    energy_threshold: float = Field(
        50.0,
        description="Minimum energy (in TWh) for an energy balance component to be included in energy plots.",
    )
    balance_timeseries: _BalanceTimeseriesConfig = Field(
        default_factory=_BalanceTimeseriesConfig,
        description="Settings for balance timeseries plots (see `plot_balance_timeseries`).",
    )
    interactive_bus_balance: _InteractiveBusBalanceConfig = Field(
        default_factory=_InteractiveBusBalanceConfig,
        description="Settings for the interactive per-bus balance plot (see `plot_interactive_bus_balance`).",
    )
    heatmap_timeseries: _HeatmapTimeseriesConfig = Field(
        default_factory=_HeatmapTimeseriesConfig,
        description="Carrier selections for marginal price, utilisation rate and state-of-charge heatmaps.",
    )
    balance_map: _BalanceMapConfig = Field(
        default_factory=_default_balance_map,
        description="Settings for static balance maps (see `plot_balance_map`).",
    )
    balance_map_interactive: _BalanceMapInteractiveConfig = Field(
        default_factory=_default_balance_map_interactive,
        description="Settings for interactive balance maps (see `plot_balance_map_interactive`).",
    )
    nice_names: dict[str, str] = Field(
        default_factory=lambda: dict(_NICE_NAMES),
        description="Human-readable display names used to relabel technology/carrier identifiers in plots and summaries.",
    )
    tech_colors: dict[str, str] = Field(
        default_factory=lambda: dict(_TECH_COLORS),
        description="Colors used to represent technologies and carriers consistently across all plots.",
    )
