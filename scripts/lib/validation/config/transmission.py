# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Transmission candidate configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#transmission
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _DynamicLineRatingConfig(ConfigModel):
    """Configuration for `lines.dynamic_line_rating` settings."""

    activate: bool = Field(
        False,
        description="Whether to take dynamic line rating into account.",
    )
    cutout: str | list[str] = Field(
        "default",
        description="Specifies the weather data cutout file(s) to use.",
    )
    correction_factor: float = Field(
        0.95,
        description="Factor to compensate for overestimation of wind speeds in hourly averaged wind data.",
    )
    max_voltage_difference: float | Literal[False] = Field(
        False,
        description="Maximum voltage angle difference in degrees or 'false' to disable.",
    )
    max_line_rating: float | Literal[False] = Field(
        False,
        description="Maximum line rating relative to nominal capacity without DLR, e.g. 1.3 or 'false' to disable.",
    )


class LinesConfig(BaseModel):
    """Configuration for `lines` settings."""

    types: dict[float, str] = Field(
        default_factory=lambda: {
            63.0: "94-AL1/15-ST1A 20.0",
            66.0: "94-AL1/15-ST1A 20.0",
            90.0: "184-AL1/30-ST1A 110.0",
            110.0: "184-AL1/30-ST1A 110.0",
            132.0: "243-AL1/39-ST1A 110.0",
            150.0: "243-AL1/39-ST1A 110.0",
            220.0: "Al/St 240/40 2-bundle 220.0",
            300.0: "Al/St 240/40 3-bundle 300.0",
            330.0: "Al/St 240/40 3-bundle 300.0",
            380.0: "Al/St 240/40 4-bundle 380.0",
            400.0: "Al/St 240/40 4-bundle 380.0",
            500.0: "Al/St 240/40 4-bundle 380.0",
            750.0: "Al/St 560/50 4-bundle 750.0",
        },
        description="Specifies line types to assume for the different voltage levels of the ENTSO-E grid extraction. Should normally handle voltage levels 220, 300, and 380 kV.",
    )
    s_max_pu: float = Field(
        0.7,
        description="Correction factor for line capacities (`s_nom`) to approximate N-1 security and reserve capacity for reactive power flows.",
    )
    s_nom_max: float = Field(
        float("inf"),
        description="Global upper limit for the maximum capacity of each extendable line (MW).",
    )
    max_extension: float = Field(
        20000,
        description="Upper limit for the extended capacity of each extendable line (MW).",
    )
    length_factor: float = Field(
        1.25,
        description="Correction factor to account for the fact that buses are *not* connected by lines through air-line distance.",
    )
    reconnect_crimea: bool = Field(
        True,
        description="Whether to reconnect Crimea to the Ukrainian grid.",
    )
    under_construction: Literal["zero", "remove", "keep"] = Field(
        "keep",
        description="Specifies how to handle lines which are currently under construction.",
    )
    dynamic_line_rating: _DynamicLineRatingConfig = Field(
        default_factory=_DynamicLineRatingConfig,
        description="Configuration for dynamic line rating.",
    )


class LinksConfig(ConfigModel):
    """Configuration for `links` settings."""

    p_max_pu: float = Field(
        1.0,
        description="Correction factor for link capacities `p_nom`.",
    )
    p_min_pu: float = Field(
        -1.0,
        description="Correction factor for link capacities `p_nom`.",
    )
    p_nom_max: float = Field(
        float("inf"),
        description="Global upper limit for the maximum capacity of each extendable DC link (MW).",
    )
    max_extension: float = Field(
        30000,
        description="Upper limit for the extended capacity of each extendable DC link (MW).",
    )
    length_factor: float = Field(
        1.25,
        description="Correction factor to account for the fact that buses are *not* connected by links through air-line distance.",
    )
    under_construction: Literal["zero", "remove", "keep"] = Field(
        "keep",
        description="Specifies how to handle lines which are currently under construction.",
    )


class TransformersConfig(ConfigModel):
    """Configuration for `transformers` settings."""

    x: float = Field(
        0.1,
        description="Series reactance in per unit (p.u.), using `s_nom` as base power of the transformer. Overwritten if `type` is specified.",
    )
    s_nom: float = Field(
        2000.0,
        description="Limit of apparent power which can pass through branch (MVA). Overwritten if `type` is specified.",
    )
    type: str = Field(
        "",
        description="Specifies transformer types to assume for the transformers of the ENTSO-E grid extraction.",
    )


class _IncludeConfig(ConfigModel):
    """Configuration for `electricity.projects.include` settings."""

    tyndp2020: bool = Field(
        True,
        description="Whether to integrate the TYNDP 2020 dataset.",
    )
    nep: bool = Field(
        True,
        description="Whether to integrate the German network development plan dataset.",
    )
    manual: bool = Field(
        True,
        description="Whether to integrate the manually added transmission projects. They are taken from the previously existing links_tyndp.csv file.",
    )


class ElectricityProjectsConfig(BaseModel):
    """Configuration for `electricity.projects` settings."""

    enable: bool = Field(
        True,
        description="Whether to integrate this transmission projects or not.",
    )
    include: _IncludeConfig = Field(
        default_factory=_IncludeConfig,
        description="Name of the transmission projects. They should be unique and have to be provided in the `data/transmission_projects` folder.",
    )
    skip: list[str] = Field(
        default_factory=lambda: ["upgraded_lines", "upgraded_links"],
        description="Type of lines to skip from all transmission projects. Possible values are: `upgraded_lines`, `upgraded_links`, `new_lines`, `new_links`.",
    )
    status: list[str] | dict[str, list[str]] = Field(
        default_factory=lambda: ["under_construction", "in_permitting", "confirmed"],
        description="Status to include into the model as list or as dict with name of project and status to include. Possible values for status are `under_construction`, `in_permitting`, `confirmed`, `planned_not_yet_permitted`, `under_consideration`.",
    )
    new_link_capacity: Literal["zero", "keep"] = Field(
        "zero",
        description="Whether to set the new link capacity to the provided capacity or set it to zero.",
    )


class _GabrielFilterConfig(BaseModel):
    """Configuration for `transmission.<carrier>.gabriel_filter` settings."""

    enable: bool = Field(
        True,
        description="Whether to filter Delaunay edges to Gabriel edges before min-degree backfilling.",
    )
    min_degree: int = Field(
        1,
        ge=0,
        description="Minimum node degree target applied after Gabriel filtering.",
    )


class _TransmissionCarrierConfigGeneral(BaseModel):
    """Configuration for a single transmission carrier."""

    enable: bool = Field(
        True,
        description="Enable transmission candidate generation for this carrier.",
    )
    gabriel_filter: _GabrielFilterConfig = Field(
        default_factory=_GabrielFilterConfig,
        description="Gabriel filter configuration.",
    )
    max_offshore_haversine_distance: float = Field(
        float("inf"),
        gt=0,
        description="Maximum haversine distance in km for offshore transmission candidates. Candidate edges a total offshore length exceeding this threshold are excluded. Defaults to infinity (no limit).",
    )
    length_factor: float = Field(
        1.25,
        gt=0,
        description="Multiplier applied to geometric corridor lengths.",
    )
    cost_factor: float = Field(
        1,
        gt=0,
        description="Multiplier applied to the capital cost of transmission infrastructure.",
    )


class _TransmissionCarrierConfigElectricity(BaseModel):
    """Configuration for electricity transmission grid."""

    enable: bool = Field(
        True,
        description="Switch for enabling/disabling the electricity transmission grid.",
    )
    base_network: Literal["entsoegridkit", "osm", "tyndp"] = Field(
        "osm",
        description="Specify the underlying base network, i.e. GridKit (based on ENTSO-E web map extract), OpenStreetMap (OSM), or TYNDP.",
    )
    transmission_limit: str = Field(
        "vopt",
        description="Limit on transmission expansion. The first part can be `v` (for setting a limit on line volume) or `c` (for setting a limit on line cost). The second part can be `opt` or a float bigger than one (e.g. 1.25). If `opt` is chosen line expansion is optimised according to its capital cost (where the choice `v` only considers overhead costs for HVDC transmission lines, while `c` uses more accurate costs distinguishing between overhead and underwater sections and including inverter pairs). The setting `v1.25` will limit the total volume of line expansion to 25% of currently installed capacities weighted by individual line lengths. The setting `c1.25` will allow to build a transmission network that costs no more than 25 % more than the current system.",
    )
    lines: LinesConfig = Field(
        default_factory=LinesConfig,
        description="Transmission lines configuration.",
    )
    links: LinksConfig = Field(
        default_factory=LinksConfig,
        description="HVDC links configuration.",
    )
    transformers: TransformersConfig = Field(
        default_factory=TransformersConfig,
        description="Transformers configuration.",
    )
    projects: ElectricityProjectsConfig = Field(
        default_factory=ElectricityProjectsConfig,
        description="Electricity transmission projects configuration.",
    )


class _TransmissionCarrierConfigGas(BaseModel):
    enable: bool = Field(
        True,
        description="Add existing natural gas infrastructure, incl. LNG terminals, production and entry-points. The existing gas network is added with a lossless transport model. A length-weighted `k-edge augmentation algorithm <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation.html#networkx.algorithms.connectivity.edge_augmentation.k_edge_augmentation>`_ can be run to add new candidate gas pipelines such that all regions of the model can be connected to the gas network. When activated, all the gas demands are regionally disaggregated as well.",
    )


class TransmissionConfig(BaseModel):
    """Configuration for `transmission` settings."""

    electricity: _TransmissionCarrierConfigElectricity = Field(
        default_factory=_TransmissionCarrierConfigElectricity,
        description="Configuration for electricity transmission grid.",
    )
    hydrogen: _TransmissionCarrierConfigGeneral = Field(
        default_factory=_TransmissionCarrierConfigGeneral,
        description="Configuration for hydrogen transmission candidates.",
    )
    carbon_dioxide: _TransmissionCarrierConfigGeneral = Field(
        default_factory=_TransmissionCarrierConfigGeneral,
        description="Configuration for carbon dioxide transmission candidates.",
    )
    gas: _TransmissionCarrierConfigGas = Field(
        default_factory=_TransmissionCarrierConfigGas,
        description="Configuration for methane gas transmission candidates.",
    )
