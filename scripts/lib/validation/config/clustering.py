# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Clustering configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#clustering
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _AdministrativeConfig(ConfigModel):
    """Configuration for `clustering.administrative` settings."""

    level: Literal[0, 1, 2, 3, "bz"] = Field(
        1,
        description="Level of administrative regions to cluster the network. 0: Country level, 1: NUTS1 level, 2: NUTS2 level, 3: NUTS3 level, 'bz': Bidding zones. Only applies when mode is set to `administrative`. Note that non-NUTS countries 'BA', 'MD', 'UA', and 'XK' can only be clustered to level 0 and 1.",
    )
    countries: dict[str, int] = Field(
        default_factory=dict,
        description="Optionally include dictionary of individual country codes and their individual NUTS levels. Overwrites country-specific `level`. For example: `{'DE': 1, 'FR': 2}`. Only applies when mode is set to `administrative`.",
    )


class _BuildBiddingZonesConfig(BaseModel):
    """Configuration for `clustering.build_bidding_zones` settings."""

    remove_islands: bool = Field(
        False,
        description="Exclude from the shape file the Balearic Islands, Bornholm, the Canary Islands, the Orkney Islands, the Shetland Islands, the Azores Islands and Madeira.",
    )
    aggregate_to_tyndp: bool = Field(
        False,
        description="Adjust the shape file to the TYNDP topology. Aggregate the Southern Norwegian bidding zones and extract Crete as a separate zone from the Greek shape.",
    )


class _SimplifyNetworkConfig(BaseModel):
    """Configuration for `clustering.simplify_network` settings."""

    to_substations: bool = Field(
        False,
        description="Aggregates all nodes without power injection (positive or negative, i.e. demand or generation) to electrically closest ones.",
    )
    exclude_carriers: list[str] = Field(
        default_factory=list,
        description="List of carriers which will not be aggregated. If empty, all carriers will be aggregated.",
    )
    remove_stubs: bool = Field(
        True,
        description="Controls whether radial parts of the network should be recursively aggregated. Defaults to true.",
    )
    remove_stubs_across_borders: bool = Field(
        False,
        description="Controls whether radial parts of the network should be recursively aggregated across borders. Defaults to true.",
    )


class _ClusterNetworkConfig(BaseModel):
    """Configuration for `clustering.cluster_network` settings."""

    algorithm: Literal["kmeans", "hac"] = Field(
        "kmeans",
        description="Clustering algorithm to use.",
    )
    hac_features: list[str] = Field(
        default_factory=lambda: ["wnd100m", "influx_direct"],
        description="List of meteorological variables contained in the weather data cutout that should be considered for hierarchical clustering.",
    )


class _AggregationStrategiesConfig(BaseModel):
    """Configuration for `clustering.aggregation_strategies` settings."""

    generators: dict[str, str] = Field(
        default_factory=lambda: {
            "committable": "any",
            "ramp_limit_up": "max",
            "ramp_limit_down": "max",
        },
        description="Aggregates the component according to the given strategy. For example, if sum, then all values within each cluster are summed to represent the new generator.",
    )
    buses: dict[str, str] = Field(
        default_factory=dict,
        description="Aggregates the component according to the given strategy. For example, if sum, then all values within each cluster are summed to represent the new bus.",
    )


class _TemporalConfig(BaseModel):
    """Configuration for `clustering.temporal` settings."""

    resolution_elec: bool | str = Field(
        False,
        description="Resample the time-resolution by averaging over every `n` snapshots in `prepare_network`. **Warning:** This option should currently only be used with electricity-only networks, not for sector-coupled networks.",
    )
    resolution_sector: bool | str = Field(
        False,
        description="Resample the time-resolution by averaging over every `n` snapshots in `prepare_sector_network`.",
    )


class ClusteringConfig(BaseModel):
    """Configuration for `clustering` settings."""

    mode: Literal["busmap", "custom_busmap", "administrative", "custom_busshapes"] = (
        Field(
            "busmap",
            description="'busmap': Default. 'custom_busmap': Enable the use of custom busmaps in rule `cluster_network`. If activated the rule looks for provided busmaps at ``data/busmaps/base_s_{clusters}_{base_network}.csv`` which should have the same format as ``resources/busmap_base_s_{clusters}.csv``, i.e. the index should contain the buses of ``networks/base_s.nc``. {base_network} is the name of the selected base_network in electricity, e.g. ``gridkit``, ``osm-prebuilt``, or ``osm-raw``. 'administrative': Clusters and indexes the network based on the administrative regions of the countries based on ``nuts3_shapes.geojson`` (level: 1, 2, 3, bz). To activate this, additionally set the ``clusters`` wildcard in ``scenario`` to 'adm'. 'custom_busshapes': Enable the use of custom shapes in rule `cluster_network`. If activated the rule looks for provided busshapes at ``data/busshapes/base_s_{clusters}_{base_network}.geojson``.",
        )
    )
    administrative: _AdministrativeConfig = Field(
        default_factory=_AdministrativeConfig,
        description="Administrative clustering settings.",
    )
    focus_weights: bool | dict[str, float] = Field(
        False,
        description="Optionally specify the focus weights for the clustering of countries. For instance: `DE: 0.8` will distribute 80% of all nodes to Germany and 20% to the rest of the countries. Only applies when mode is set to `busmap`.",
    )
    copperplate_regions: list[list[str]] = Field(
        default_factory=list,
        description="Optionally specify the regions to copperplate as a list of groups. Each group is a list of region codes that will be connected with infinite capacity lines.",
    )
    build_bidding_zones: _BuildBiddingZonesConfig = Field(
        default_factory=_BuildBiddingZonesConfig,
        description="Build bidding zones configuration.",
    )
    simplify_network: _SimplifyNetworkConfig = Field(
        default_factory=_SimplifyNetworkConfig,
        description="Network simplification settings.",
    )
    cluster_network: _ClusterNetworkConfig = Field(
        default_factory=_ClusterNetworkConfig,
        description="Network clustering algorithm settings.",
    )
    exclude_carriers: list[str] = Field(
        default_factory=list,
        description="List of carriers which will not be aggregated. If empty, all carriers will be aggregated.",
    )
    consider_efficiency_classes: bool = Field(
        False,
        description="Aggregated each carriers into the top 10-quantile (high), the bottom 90-quantile (low), and everything in between (medium).",
    )
    aggregation_strategies: _AggregationStrategiesConfig = Field(
        default_factory=_AggregationStrategiesConfig,
        description="Aggregation strategies for different components.",
    )
    temporal: _TemporalConfig = Field(
        default_factory=_TemporalConfig,
        description="Options for temporal resolution.",
    )
