# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build coefficient of performance (COP) time series for air- or ground-sourced
heat pumps.

For individual (decentral) heat pumps, the COP is approximated as a quatratic function of the temperature difference between source and sink, based on Staffell et al. 2012.
For district (central) heating, the COP is approximated based on Jensen et al. 2018 and parameters from Pieper et al. 2020.

This rule is executed in ``build_sector.smk``.

Relevant Settings
-----------------

.. code:: yaml
    heat_pump_sink_T:


Inputs:
-------
- ``resources/<run_name>/temp_soil_total_elec_s<simpl>_<clusters>.nc``: Soil temperature (total) time series.
- ``resources/<run_name>/temp_air_total_elec_s<simpl>_<clusters>.nc``: Ambient air temperature (total) time series.

Outputs:
--------
- ``resources/cop_air_decentral_heating_elec_s<simpl>_<clusters>.nc``: COP (air-sourced) time series (decentral heating).
- ``resources/cop_soil_decentral_heating_elec_s<simpl>_<clusters>.nc``: COP (ground-sourced) time series (decentral heating).
- ``resources/cop_air_central_heating_elec_s<simpl>_<clusters>.nc``: COP (air-sourced) time series (central heating).
- ``resources/cop_soil_central_heating_elec_s<simpl>_<clusters>.nc``: COP (ground-sourced) time series (central heating).


References
----------
[1] Staffell et al., Energy & Environmental Science 11 (2012): A review of domestic heat pumps, https://doi.org/10.1039/C2EE22653G.
[2] Jensen et al., Proceedings of the13th IIR-Gustav Lorentzen Conference on Natural Refrigerants (2018): Heat pump COP, part 2: Generalized COP estimation of heat pump processes, https://doi.org/10.18462/iir.gl.2018.1386
[3] Pieper et al., Energy 205 (2020): Comparison of COP estimation methods for large-scale heat pumps used in energy planning, https://doi.org/10.1016/j.energy.2020.117994
"""

from enum import Enum
from typing import Union

import numpy as np
import xarray as xr
from _helpers import set_scenario_config


def coefficient_of_performance_individual_heating(delta_T, source="air"):
    if source == "air":
        return 6.81 - 0.121 * delta_T + 0.000630 * delta_T**2
    elif source == "soil":
        return 8.77 - 0.150 * delta_T + 0.000734 * delta_T**2
    else:
        raise NotImplementedError("'source' must be one of  ['air', 'soil']")


def celsius_to_kelvin(
    t_celsius: Union[float, xr.DataArray, np.array]
) -> Union[float, xr.DataArray, np.array]:
    if (np.asarray(t_celsius) > 200).any():
        raise ValueError("t_celsius > 200. Are you sure you are using the right units?")
    return t_celsius + 273.15


def logarithmic_mean(
    t_hot: Union[float, xr.DataArray, np.ndarray],
    t_cold: Union[float, xr.DataArray, np.ndarray],
) -> Union[float, xr.DataArray, np.ndarray]:
    if (np.asarray(t_hot <= t_cold)).any():
        raise ValueError("t_hot must be greater than t_cold")
    return (t_hot - t_cold) / np.log(t_hot / t_cold)


class CopDistrictHeating:

    def __init__(
        self,
        forward_temperature_celsius: Union[xr.DataArray, np.array],
        source_inlet_temperature_celsius: Union[xr.DataArray, np.array],
        return_temperature_celsius: Union[xr.DataArray, np.array],
        source_outlet_temperature_celsius: Union[xr.DataArray, np.array],
        delta_t_pinch_point: float = 5,
        isentropic_compressor_efficiency: float = 0.8,
        heat_loss: float = 0.0,
    ) -> None:
        """
        Initialize the COPProfileBuilder object.

        Parameters:
        ----------
        forward_temperature_celsius : Union[xr.DataArray, np.array]
            The forward temperature in Celsius.
        return_temperature_celsius : Union[xr.DataArray, np.array]
            The return temperature in Celsius.
        source_inlet_temperature_celsius : Union[xr.DataArray, np.array]
            The source inlet temperature in Celsius.
        source_outlet_temperature_celsius : Union[xr.DataArray, np.array]
            The source outlet temperature in Celsius.
        delta_t_pinch_point : float, optional
            The pinch point temperature difference, by default 5.
        isentropic_compressor_efficiency : float, optional
            The isentropic compressor efficiency, by default 0.8.
        heat_loss : float, optional
            The heat loss, by default 0.0.
        """
        self.t_source_in = celsius_to_kelvin(source_inlet_temperature_celsius)
        self.t_sink_out = celsius_to_kelvin(forward_temperature_celsius)

        self.t_sink_in = celsius_to_kelvin(return_temperature_celsius)
        self.t_source_out = celsius_to_kelvin(source_outlet_temperature_celsius)

        self.isentropic_efficiency_compressor = isentropic_compressor_efficiency
        self.heat_loss = heat_loss
        self.delta_t_pinch = delta_t_pinch_point

    def cop(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the coefficient of performance (COP) for the system.

        Returns:
            Union[xr.DataArray, np.array]: The calculated COP values.
        """
        return (
            self.ideal_lorenz_cop
            * (
                (
                    1
                    + (self.delta_t_refrigerant_sink + self.delta_t_pinch)
                    / self.t_sink_mean
                )
                / (
                    1
                    + (
                        self.delta_t_refrigerant_sink
                        + self.delta_t_refrigerant_source
                        + 2 * self.delta_t_pinch
                    )
                    / self.delta_t_lift
                )
            )
            * self.isentropic_efficiency_compressor
            * (1 - self.ratio_evaporation_compression_work)
            + 1
            - self.isentropic_efficiency_compressor
            - self.heat_loss
        )

    @property
    def t_sink_mean(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the logarithmic mean temperature difference between the cold
        and hot sinks.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The mean temperature difference.
        """
        return logarithmic_mean(t_cold=self.t_sink_in, t_hot=self.t_sink_out)

    @property
    def t_source_mean(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the logarithmic mean temperature of the heat source.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The mean temperature of the heat source.
        """
        return logarithmic_mean(t_hot=self.t_source_in, t_cold=self.t_source_out)

    @property
    def delta_t_lift(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the temperature lift as the difference between the
        logarithmic sink and source temperatures.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The temperature difference between the sink and source.
        """
        return self.t_sink_mean - self.t_source_mean

    @property
    def ideal_lorenz_cop(self) -> Union[xr.DataArray, np.array]:
        """
        Ideal Lorenz coefficient of performance (COP).

        The ideal Lorenz COP is calculated as the ratio of the mean sink temperature
        to the lift temperature difference.

        Returns
        -------
        np.array
            The ideal Lorenz COP.
        """
        return self.t_sink_mean / self.delta_t_lift

    @property
    def delta_t_refrigerant_source(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the temperature difference between the refrigerant source
        inlet and outlet.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The temperature difference between the refrigerant source inlet and outlet.
        """
        return self._approximate_delta_t_refrigerant_source(
            delta_t_source=self.t_source_in - self.t_source_out
        )

    @property
    def delta_t_refrigerant_sink(self) -> Union[xr.DataArray, np.array]:
        """
        Temperature difference between the refrigerant and the sink based on
        approximation.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The temperature difference between the refrigerant and the sink.
        """
        return self._approximate_delta_t_refrigerant_sink()

    @property
    def ratio_evaporation_compression_work(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the ratio of evaporation to compression work based on
        approximation.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The calculated ratio of evaporation to compression work.
        """
        return self._ratio_evaporation_compression_work_approximation()

    @property
    def delta_t_sink(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the temperature difference at the sink.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The temperature difference at the sink.
        """
        return self.t_sink_out - self.t_sink_in

    def _approximate_delta_t_refrigerant_source(
        self, delta_t_source: Union[xr.DataArray, np.array]
    ) -> Union[xr.DataArray, np.array]:
        """
        Approximates the temperature difference between the refrigerant and the
        source.

        Parameters
        ----------
        delta_t_source : Union[xr.DataArray, np.array]
            The temperature difference for the refrigerant source.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The approximate temperature difference for the refrigerant source.
        """
        return delta_t_source / 2

    def _approximate_delta_t_refrigerant_sink(
        self,
        refrigerant: str = "ammonia",
        a: float = {"ammonia": 0.2, "isobutane": -0.0011},
        b: float = {"ammonia": 0.2, "isobutane": 0.3},
        c: float = {"ammonia": 0.016, "isobutane": 2.4},
    ) -> Union[xr.DataArray, np.array]:
        """
        Approximates the temperature difference at the refrigerant sink.

        Parameters:
        ----------
        refrigerant : str, optional
            The refrigerant used in the system. Either 'isobutane' or 'ammonia. Default is 'ammonia'.
        a : float, optional
            Coefficient for the temperature difference between the sink and source, default is 0.2.
        b : float, optional
            Coefficient for the temperature difference at the sink, default is 0.2.
        c : float, optional
            Constant term, default is 0.016.

        Returns:
        -------
        Union[xr.DataArray, np.array]
            The approximate temperature difference at the refrigerant sink.

        Notes:
        ------
        This function assumes ammonia as the refrigerant.

        The approximate temperature difference at the refrigerant sink is calculated using the following formula:
        a * (t_sink_out - t_source_out + 2 * delta_t_pinch) + b * delta_t_sink + c
        """
        if refrigerant not in a.keys():
            raise ValueError(
                f"Invalid refrigerant '{refrigerant}'. Must be one of {a.keys()}"
            )
        return (
            a[refrigerant]
            * (self.t_sink_out - self.t_source_out + 2 * self.delta_t_pinch)
            + b[refrigerant] * self.delta_t_sink
            + c[refrigerant]
        )

    def _ratio_evaporation_compression_work_approximation(
        self,
        refrigerant: str = "ammonia",
        a: float = {"ammonia": 0.0014, "isobutane": 0.0035},
        b: float = {"ammonia": -0.0015, "isobutane": -0.0033},
        c: float = {"ammonia": 0.039, "isobutane": 0.053},
    ) -> Union[xr.DataArray, np.array]:
        """
        Calculate the ratio of evaporation to compression work approximation.

        Parameters:
        ----------
        refrigerant : str, optional
            The refrigerant used in the system. Either 'isobutane' or 'ammonia. Default is 'ammonia'.
        a : float, optional
            Coefficient 'a' in the approximation equation. Default is 0.0014.
        b : float, optional
            Coefficient 'b' in the approximation equation. Default is -0.0015.
        c : float, optional
            Coefficient 'c' in the approximation equation. Default is 0.039.

        Returns:
        -------
        Union[xr.DataArray, np.array]
            The calculated ratio of evaporation to compression work.

        Notes:
        ------
        This function assumes ammonia as the refrigerant.

        The approximation equation used is:
        ratio = a * (t_sink_out - t_source_out + 2 * delta_t_pinch) + b * delta_t_sink + c
        """
        if refrigerant not in a.keys():
            raise ValueError(
                f"Invalid refrigerant '{refrigerant}'. Must be one of {a.keys()}"
            )
        return (
            a[refrigerant]
            * (self.t_sink_out - self.t_source_out + 2 * self.delta_t_pinch)
            + b[refrigerant] * self.delta_t_sink
            + c[refrigerant]
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            simpl="",
            clusters=48,
        )

    set_scenario_config(snakemake)

    for source in ["air", "soil"]:
        source_T = xr.open_dataarray(snakemake.input[f"temp_{source}_total"])

        delta_T = snakemake.params.heat_pump_sink_T_decentral_heating - source_T

        cop_individual_heating = coefficient_of_performance_individual_heating(
            delta_T, source
        )
        cop_individual_heating.to_netcdf(
            snakemake.output[f"cop_{source}_decentral_heating"]
        )

        cop_district_heating = CopDistrictHeating(
            forward_temperature_celsius=snakemake.params.forward_temperature_district_heating,
            return_temperature_celsius=snakemake.params.return_temperature_district_heating,
            source_inlet_temperature_celsius=source_T,
            source_outlet_temperature_celsius=source_T
            - snakemake.params.heat_source_cooling_district_heating,
        ).cop()

        cop_district_heating.to_netcdf(
            snakemake.output[f"cop_{source}_central_heating"]
        )
