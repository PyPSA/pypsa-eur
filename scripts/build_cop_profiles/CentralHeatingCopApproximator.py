# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


from typing import Union

import numpy as np
import xarray as xr
from BaseCopApproximator import BaseCopApproximator


class CentralHeatingCopApproximator(BaseCopApproximator):
    """
    Approximate the coefficient of performance (COP) for a heat pump in a
    central heating system (district heating).

    Uses an approximation method proposed by Jensen et al. (2018) and
    default parameters from Pieper et al. (2020). The method is based on
    a thermodynamic heat pump model with some hard-to-know parameters
    being approximated.

    Attributes
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

    Methods
    -------
    __init__(
        forward_temperature_celsius: Union[xr.DataArray, np.array],
        source_inlet_temperature_celsius: Union[xr.DataArray, np.array],
        return_temperature_celsius: Union[xr.DataArray, np.array],
        source_outlet_temperature_celsius: Union[xr.DataArray, np.array],
        delta_t_pinch_point: float = 5,
        isentropic_compressor_efficiency: float = 0.8,
        heat_loss: float = 0.0,
    ) -> None:
        Initializes the CentralHeatingCopApproximator object.

    approximate_cop(self) -> Union[xr.DataArray, np.array]:
        Calculate the coefficient of performance (COP) for the system.

    _approximate_delta_t_refrigerant_source(
        self, delta_t_source: Union[xr.DataArray, np.array]
    ) -> Union[xr.DataArray, np.array]:
        Approximates the temperature difference between the refrigerant and the source.

    _approximate_delta_t_refrigerant_sink(
        self,
        refrigerant: str = "ammonia",
        a: float = {"ammonia": 0.2, "isobutane": -0.0011},
        b: float = {"ammonia": 0.2, "isobutane": 0.3},
        c: float = {"ammonia": 0.016, "isobutane": 2.4},
    ) -> Union[xr.DataArray, np.array]:
        Approximates the temperature difference between the refrigerant and heat sink.

    _ratio_evaporation_compression_work_approximation(
        self,
        refrigerant: str = "ammonia",
        a: float = {"ammonia": 0.0014, "isobutane": 0.0035},
    ) -> Union[xr.DataArray, np.array]:
        Calculate the ratio of evaporation to compression work based on approximation.

    _approximate_delta_t_refrigerant_sink(
        self,
        refrigerant: str = "ammonia",
        a: float = {"ammonia": 0.2, "isobutane": -0.0011},
        b: float = {"ammonia": 0.2, "isobutane": 0.3},
        c: float = {"ammonia": 0.016, "isobutane": 2.4},
    ) -> Union[xr.DataArray, np.array]:
        Approximates the temperature difference between the refrigerant and heat sink.

    _ratio_evaporation_compression_work_approximation(
        self,
        refrigerant: str = "ammonia",
        a: float = {"ammonia": 0.0014, "isobutane": 0.0035},
    ) -> Union[xr.DataArray, np.array]:
        Calculate the ratio of evaporation to compression work based on approximation.
    """

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
        Initializes the CentralHeatingCopApproximator object.

        Parameters
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
        self.t_source_in_kelvin = BaseCopApproximator.celsius_to_kelvin(
            source_inlet_temperature_celsius
        )
        self.t_sink_out_kelvin = BaseCopApproximator.celsius_to_kelvin(
            forward_temperature_celsius
        )

        self.t_sink_in_kelvin = BaseCopApproximator.celsius_to_kelvin(
            return_temperature_celsius
        )
        self.t_source_out = BaseCopApproximator.celsius_to_kelvin(
            source_outlet_temperature_celsius
        )

        self.isentropic_efficiency_compressor_kelvin = isentropic_compressor_efficiency
        self.heat_loss = heat_loss
        self.delta_t_pinch = delta_t_pinch_point

    def approximate_cop(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the coefficient of performance (COP) for the system.

        Notes
        -----
        Returns 0 where the source inlet temperature is greater than the sink outlet temperature.

        Returns
        -------
            Union[xr.DataArray, np.array]: The calculated COP values.
        """
        return xr.where(
            self.t_source_in_kelvin > self.t_sink_out_kelvin,
            0,
            self.ideal_lorenz_cop
            * (
                (
                    1
                    + (self.delta_t_refrigerant_sink + self.delta_t_pinch)
                    / self.t_sink_mean_kelvin
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
            * self.isentropic_efficiency_compressor_kelvin
            * (1 - self.ratio_evaporation_compression_work)
            + 1
            - self.isentropic_efficiency_compressor_kelvin
            - self.heat_loss,
        )

    @property
    def t_sink_mean_kelvin(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the logarithmic mean temperature difference between the cold
        and hot sinks.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The mean temperature difference.
        """
        return BaseCopApproximator.logarithmic_mean(
            t_cold=self.t_sink_in_kelvin, t_hot=self.t_sink_out_kelvin
        )

    @property
    def t_source_mean_kelvin(self) -> Union[xr.DataArray, np.array]:
        """
        Calculate the logarithmic mean temperature of the heat source.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The mean temperature of the heat source.
        """
        return BaseCopApproximator.logarithmic_mean(
            t_hot=self.t_source_in_kelvin, t_cold=self.t_source_out
        )

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
        return self.t_sink_mean_kelvin - self.t_source_mean_kelvin

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
        return self.t_sink_mean_kelvin / self.delta_t_lift

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
            delta_t_source=self.t_source_in_kelvin - self.t_source_out
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
        return self.t_sink_out_kelvin - self.t_sink_in_kelvin

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
            The approximate temperature difference between the refrigerant and heat source.
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
        Approximates the temperature difference between the refrigerant and
        heat sink.

        Parameters
        ----------
        refrigerant : str, optional
            The refrigerant used in the system. Either 'isobutane' or 'ammonia. Default is 'ammonia'.
        a : float, optional
            Coefficient for the temperature difference between the sink and source, default is 0.2.
        b : float, optional
            Coefficient for the temperature difference at the sink, default is 0.2.
        c : float, optional
            Constant term, default is 0.016.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The approximate temperature difference between the refrigerant and heat sink.

        Notes
        -----
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
            * (self.t_sink_out_kelvin - self.t_source_out + 2 * self.delta_t_pinch)
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

        Parameters
        ----------
        refrigerant : str, optional
            The refrigerant used in the system. Either 'isobutane' or 'ammonia. Default is 'ammonia'.
        a : float, optional
            Coefficient 'a' in the approximation equation. Default is 0.0014.
        b : float, optional
            Coefficient 'b' in the approximation equation. Default is -0.0015.
        c : float, optional
            Coefficient 'c' in the approximation equation. Default is 0.039.

        Returns
        -------
        Union[xr.DataArray, np.array]
            The approximated ratio of evaporation to compression work.

        Notes
        -----
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
            * (self.t_sink_out_kelvin - self.t_source_out + 2 * self.delta_t_pinch)
            + b[refrigerant] * self.delta_t_sink
            + c[refrigerant]
        )
