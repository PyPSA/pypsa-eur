# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculate the Aquifer Thermal Energy Storage (ATES) potentials for each region based on suitable aquifers, district heating areas, and temperature profiles.

The script filters aquifers based on their geological suitability (currently only highly productive porous aquifers), intersects them with future district heating areas and onshore regions, and calculates the energy storage potential. The difference between average forward and return temperatures is used as temperature differentials to convert from volumetric to energetic storage potentials.

The methodology is loosely based on Jackson, Regnier, Staffell (2024).
Future district heating areas are sourced from Manz et al. (2024), based on Fallahnejad et al. (2024). Data on aquifers stems from the German Bundesanstalt für Geowissenschaften und Rohstoffe (BGR).


Relevant Settings
-----------------

.. code:: yaml
    sector:
        aquifer_thermal_energy_storage:
            aquifer_volumetric_heat_capacity:
            fraction_of_aquifer_area_available:
            effective_screen_length:
            suitable_aquifer_types:
            dh_area_buffer:

Inputs
------
- `resources/<run_name>/regions_onshore.geojson`: Shapes of onshore regions
- `resources/<run_name>/aquifer_shapes.shp`: Shapes of aquifers
- `resources/<run_name>/dh_areas.geojson`: Shapes of district heating areas
- `resources/<run_name>/central_heating_forward_temperature_profiles.nc`: Forward temperature profiles
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`: Return temperature profiles

Outputs
-------
- `resources/<run_name>/ates_potentials.csv`: ATES potentials per region in MWh

References
----------
- Jackson, Regnier, Staffell 2024: "Aquifer Thermal Energy Storage for low carbon heating and cooling in the United Kingdom: Current status and future prospects", Applied Energy, vol. 376, no. 124096, https://doi.org/10.1016/j.apenergy.2024.124096
- Manz et al. 2024: "Spatial analysis of renewable and excess heat potentials for climate-neutral district heating in Europe", Renewable Energy, vol. 224, no. 120111, https://doi.org/10.1016/j.renene.2024.120111
- Fallahnejad et al. 2024: "District heating potential in the EU-27: Evaluating the impacts of heat demand reduction and market share growth", Applied Energy, vol. 353, no. 122154, https://https://doi.org/10.1016/j.apenergy.2023.122154
- BGR: IHME1500 - Internationale Hydrogeologische Karte von Europa 1:1.500.000 (https://www.bgr.bund.de/DE/Themen/Wasser/Projekte/laufend/Beratung/Ihme1500/ihme1500_projektbeschr.html?nn=1546102)
"""

import logging
import sys
import os
from pathlib import Path

import geopandas as gpd
import pandas as pd
import xarray as xr
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def mwh_ates_per_m2(
    aquifer_volumetric_heat_capacity: float,
    fraction_of_aquifer_area_available: float,
    effective_screen_length: float,
    hot_well_temperature: float,
    cold_well_temperature: float,
    kwh_per_kj: float = 1 / 3600,
    mwh_per_kwh: float = 1 / 1e3,
) -> float:
    """
    Calculate the MWh potential per square meter for ATES systems.

    This is based on Jackson, Regnier, Staffell 2024 (https://doi.org/10.1016/j.apenergy.2024.124096).
    
    Parameters
    ----------
    aquifer_volumetric_heat_capacity : float
        Volumetric heat capacity of the aquifer in kJ/(m³·K)
    fraction_of_aquifer_area_available : float
        Fraction of the aquifer area that can be used for ATES (between 0 and 1)
    effective_screen_length : float
        Effective depth/thickness of the aquifer in meters
    hot_well_temperature : float
        Temperature of the hot well in degrees Celsius
    cold_well_temperature : float
        Temperature of the cold well in degrees Celsius
    kwh_per_kj : float, optional
        Conversion factor from kJ to kWh, by default 1/3600
    mwh_per_kwh : float, optional
        Conversion factor from kWh to MWh, by default 1/1000
    
    Returns
    -------
    float
        The ATES energy storage potential in MWh per square meter
        
    Raises
    ------
    Exception
        If calculation fails
    
    References
    ----------
    - Jackson, Regnier, Staffell 2024 (https://doi.org/10.1016/j.apenergy.2024.124096): Aquifer Thermal Energy Storages for low carbon heating and cooling in the United Kingdom: CUrrent status and future prospects
    """
    try:
        return (
            aquifer_volumetric_heat_capacity
            * fraction_of_aquifer_area_available
            * effective_screen_length
            * (hot_well_temperature - cold_well_temperature)
            * kwh_per_kj
            * mwh_per_kwh
        )
    except Exception as e:
        logger.error(f"Error calculating ATES potential per m2: {e}")
        raise


def suitable_aquifers(
    aquifer_shapes: gpd.GeoDataFrame,
    suitable_aquifer_types: list,
) -> gpd.GeoDataFrame:
    """
    Filter the aquifer shapes by the suitable aquifer types.
    
    Parameters
    ----------
    aquifer_shapes : geopandas.GeoDataFrame
        GeoDataFrame containing the shapes of all aquifers
    suitable_aquifer_types : list
        List of aquifer types considered suitable for ATES
    
    Returns
    -------
    geopandas.GeoDataFrame
        Filtered GeoDataFrame containing only suitable aquifers
        
    Raises
    ------
    KeyError
        If required column 'AQUIF_NAME' is not found in the input data
    Exception
        If filtering process fails
    """
    try:
        if "AQUIF_NAME" not in aquifer_shapes.columns:
            raise KeyError("Column 'AQUIF_NAME' not found in aquifer shapes dataframe")
        
        filtered = aquifer_shapes[aquifer_shapes["AQUIF_NAME"].isin(suitable_aquifer_types)]
        
        if filtered.empty:
            logger.warning("No suitable aquifers found with the specified types")
            
        return filtered
    except Exception as e:
        logger.error(f"Error filtering suitable aquifers: {e}")
        raise


def ates_potential_per_onshore_region(
    suitable_aquifers: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    dh_areas: gpd.GeoDataFrame,
    dh_area_buffer: float,
    mwh_per_m2: float,
) -> gpd.GeoDataFrame:
    """
    Calculate the ATES potential for each onshore region.
    
    This function overlays suitable aquifers with onshore regions and district heating areas
    to calculate the potential ATES capacity for each region.
    
    Parameters
    ----------
    suitable_aquifers : geopandas.GeoDataFrame
        GeoDataFrame containing filtered suitable aquifers
    regions_onshore : geopandas.GeoDataFrame
        GeoDataFrame containing the shapes of onshore regions
    dh_areas : geopandas.GeoDataFrame
        GeoDataFrame containing the shapes of district heating areas
    dh_area_buffer : float
        Buffer distance in meters to apply around district heating areas
    mwh_per_m2 : float
        ATES potential in MWh per square meter
    
    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame with onshore regions and their ATES potential in MWh
        
    Raises
    ------
    KeyError
        If required column 'name' is not found in regions_onshore dataframe
    Exception
        If calculation process fails
    """
    try:
        if suitable_aquifers.empty or regions_onshore.empty or dh_areas.empty:
            logger.warning("One or more input GeoDataFrames are empty")
            
        ret_val = regions_onshore.copy(deep=True)
        
        if "name" not in ret_val.columns:
            raise KeyError("Column 'name' not found in regions_onshore dataframe")
            
        ret_val.index = ret_val["name"]
        ret_val.drop(columns=["name"], inplace=True)

        suitable_aquifers_in_onshore_regions = gpd.overlay(
            suitable_aquifers, regions_onshore, how="intersection"
        )

        if suitable_aquifers_in_onshore_regions.empty:
            logger.warning("No suitable aquifers found in onshore regions")
            ret_val["ates_potential"] = 0
            return ret_val

        dh_areas_buffered = dh_areas.copy(deep=True)
        dh_areas_buffered["geometry"] = dh_areas_buffered.geometry.buffer(dh_area_buffer)

        try:
            aquifers_in_dh_areas = (
                gpd.overlay(
                    dh_areas_buffered, suitable_aquifers_in_onshore_regions, how="intersection"
                )
                .groupby("name")["geometry"]
                .apply(lambda x: x.area.sum())
            )
            
            # Handle regions without any ATES potential
            missing_regions = set(ret_val.index) - set(aquifers_in_dh_areas.index)
            if missing_regions:
                logger.info(f"{len(missing_regions)} regions have no ATES potential")
                
            ret_val["ates_potential"] = 0  # Default value
            ret_val.loc[aquifers_in_dh_areas.index, "ates_potential"] = aquifers_in_dh_areas * mwh_per_m2
            
        except Exception as e:
            logger.error(f"Error in overlay calculation: {e}")
            ret_val["ates_potential"] = 0
            
        return ret_val
        
    except Exception as e:
        logger.error(f"Error calculating ATES potential per onshore region: {e}")
        raise


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake("build_ates_potential", clusters="48")
        
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # get onshore regions and index them by region name
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)
    regions_onshore = regions_onshore.to_crs(epsg=3035)

    aquifer_shapes = gpd.read_file(snakemake.input.aquifer_shapes_shp).to_crs(regions_onshore.crs)

    dh_areas = gpd.read_file(snakemake.input.dh_areas).to_crs(regions_onshore.crs)

    forward_temperature_profiles = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles)
    return_temperature_profiles = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    logger.info("Calculating ATES potentials per m2")
    mwh_per_m2 = mwh_ates_per_m2(
        aquifer_volumetric_heat_capacity=snakemake.params.aquifer_volumetric_heat_capacity,
        fraction_of_aquifer_area_available=snakemake.params.fraction_of_aquifer_area_available,
        effective_screen_length=snakemake.params.effective_screen_length,
        hot_well_temperature=forward_temperature_profiles.mean(dim="time").to_pandas(),
        cold_well_temperature=return_temperature_profiles.to_pandas(),
    )

    logger.info("Filtering for suitable aquifers")
    suitable_aquifer_df = suitable_aquifers(
        aquifer_shapes=aquifer_shapes,
        suitable_aquifer_types=snakemake.params.suitable_aquifer_types,
    )

    logger.info("Calculating ATES potentials per region")
    ates_potentials = ates_potential_per_onshore_region(
        suitable_aquifers=suitable_aquifer_df,
        regions_onshore=regions_onshore,
        dh_areas=dh_areas,
        dh_area_buffer=snakemake.params.dh_area_buffer,
        mwh_per_m2=mwh_per_m2,
    )

    logger.info(f"Writing results to {snakemake.output.ates_potentials}")
    ates_potentials.to_csv(snakemake.output.ates_potentials, index_label="name")

