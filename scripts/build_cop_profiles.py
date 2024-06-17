# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build coefficient of performance (COP) time series for air- or ground-sourced
heat pumps.

The COP is approximated as a quatratic function of the temperature difference between source and
sink, based on Staffell et al. 2012.

This rule is executed in ``build_sector.smk``.

Relevant Settings
-----------------

.. code:: yaml
    heat_pump_sink_T:


Inputs:
-------
- ``resources/<run_name>/temp_soil_total_elec_s<simpl>_<clusters>.nc``: Soil temperature (total) time series.
- ``resources/<run_name>/temp_soil_rural_elec_s<simpl>_<clusters>.nc``: Soil temperature (rural) time series.
- ``resources/<run_name>/temp_soil_urban_elec_s<simpl>_<clusters>.nc``: Soil temperature (urban) time series.
- ``resources/<run_name>/temp_air_total_elec_s<simpl>_<clusters>.nc``: Ambient air temperature (total) time series.
- ``resources/<run_name>/temp_air_rural_elec_s<simpl>_<clusters>.nc``: Ambient air temperature (rural) time series.
- ``resources/<run_name>/temp_air_urban_elec_s<simpl>_<clusters>.nc``: Ambient air temperature (urban) time series.

Outputs:
--------
- ``resources/cop_soil_total_elec_s<simpl>_<clusters>.nc``: COP (ground-sourced) time series (total).
- ``resources/cop_soil_rural_elec_s<simpl>_<clusters>.nc``: COP (ground-sourced) time series (rural).
- ``resources/cop_soil_urban_elec_s<simpl>_<clusters>.nc``: COP (ground-sourced) time series (urban).
- ``resources/cop_air_total_elec_s<simpl>_<clusters>.nc``: COP (air-sourced) time series (total).
- ``resources/cop_air_rural_elec_s<simpl>_<clusters>.nc``: COP (air-sourced) time series (rural).
- ``resources/cop_air_urban_elec_s<simpl>_<clusters>.nc``: COP (air-sourced) time series (urban).


References
----------
[1] Staffell et al., Energy & Environmental Science 11 (2012): A review of domestic heat pumps, https://doi.org/10.1039/C2EE22653G.
"""

import xarray as xr
from _helpers import set_scenario_config


def coefficient_of_performance(delta_T, source="air"):
    if source == "air":
        return 6.81 - 0.121 * delta_T + 0.000630 * delta_T**2
    elif source == "soil":
        return 8.77 - 0.150 * delta_T + 0.000734 * delta_T**2
    else:
        raise NotImplementedError("'source' must be one of  ['air', 'soil']")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            simpl="",
            clusters=48,
        )

    set_scenario_config(snakemake)

    for area in ["total", "urban", "rural"]:
        for source in ["air", "soil"]:
            source_T = xr.open_dataarray(snakemake.input[f"temp_{source}_{area}"])

            delta_T = snakemake.params.heat_pump_sink_T - source_T

            cop = coefficient_of_performance(delta_T, source)

            cop.to_netcdf(snakemake.output[f"cop_{source}_{area}"])
