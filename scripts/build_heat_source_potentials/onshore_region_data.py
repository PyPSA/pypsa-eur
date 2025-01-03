# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Helper class for matching heat source potentials to onshore regions.
"""

import geopandas as gpd


class OnshoreRegionData:
    """
    This class is used to map heat potentials to onshore regions.

    Attributes
    ----------
    onshore_regions : gpd.GeoDataFrame
        GeoDataFrame containing the onshore regions
    data : gpd.GeoDataFrame
        GeoDataFrame containing the heat potentials
    scaling_factor : float
        Scaling factor for the heat potentials
    """

    def __init__(
        self,
        onshore_regions: gpd.GeoDataFrame,
        data: gpd.GeoDataFrame,
        column_name: str,
        scaling_factor: float = 1.0,
    ) -> None:
        """
        Parameters
        ----------
        onshore_regions : gpd.GeoDataFrame
            GeoDataFrame containing the onshore regions
        data : gpd.GeoDataFrame
            GeoDataFrame containing the heat potentials
        column_name : str
            Column name of the heat potential data in `data`
        scaling_factor : float, optional
            Scaling factor for the heat potentials, by default 1.0
        """

        self.onshore_regions = onshore_regions
        self.scaling_factor = scaling_factor
        self.data = data.to_crs(onshore_regions.crs)
        self._column_name = column_name

        self._mapped = False

    @property
    def data_in_regions_scaled(self) -> gpd.GeoDataFrame:
        """
        Scale the heat potentials and map them to the onshore regions.

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame containing the scaled heat potentials in the onshore regions
        """
        if self._mapped:
            return self._scaled_data_in_regions
        else:
            self._data_in_regions = self._map_to_onshore_regions()
            self._mapped = True
            return self._scaled_data_in_regions

    def _map_to_onshore_regions(self):
        """
        Map the heat potentials to the onshore regions
        """
        data_in_regions = gpd.sjoin(self.data, self.onshore_regions, how="right")

        # Initialize an empty list to store the merged GeoDataFrames
        ret_val = self.onshore_regions.copy()
        ret_val[self._column_name] = (
            data_in_regions.groupby("name")[self._column_name]
            .sum()
            .reset_index(drop=True)
        )
        ret_val = ret_val.set_index("name", drop=True).rename_axis("name")[
            self._column_name
        ]

        return ret_val

    @property
    def _scaled_data_in_regions(self):
        """
        Scale the heat potentials in the onshore regions
        """
        return self._data_in_regions * self.scaling_factor
