# -*- coding: utf-8 -*-
from typing import List

import geopandas as gpd


class OnshoreRegionData:

    def __init__(
        self,
        onshore_regions: gpd.GeoDataFrame,
        data: gpd.GeoDataFrame,
        column_name: str,
        scaling_factor: float = 1.0,
    ) -> None:

        self.onshore_regions = onshore_regions
        self.column_name = column_name
        self.scaling_factor = scaling_factor
        self.data = data.to_crs(onshore_regions.crs)

        self._mapped = False

    @property
    def data_in_regions_scaled(self) -> gpd.GeoDataFrame:
        if self._mapped:
            return self._scaled_data_in_regions
        else:
            self._data_in_regions = self._map_to_onshore_regions()
            self._mapped = True
            return self._scaled_data_in_regions

    def _map_to_onshore_regions(self):
        """
        This function maps the heat potentials to the onshore regions


        """
        data_in_regions = gpd.sjoin(self.data, self.onshore_regions, how="right")

        # Initialize an empty list to store the merged GeoDataFrames
        ret_val = self.onshore_regions.copy()
        ret_val[self.column_name] = (
            data_in_regions.groupby("name")[self.column_name]
            .sum()
            .reset_index(drop=True)
        )
        ret_val = ret_val.set_index("name", drop=True).rename_axis("name")[
            self.column_name
        ]

        return ret_val

    @property
    def _scaled_data_in_regions(self):
        # scaled_data_in_regions = self._data_in_regions.copy()
        # scaled_data_in_regions[self.column_name] *= self.scaling_factor
        return self._data_in_regions * self.scaling_factor
