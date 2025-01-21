# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class CentralHeatingTemperatureApproximator:
    """
    A class to approximate central heating temperatures based on ambient
    temperature.

    Attributes
    ----------
    ambient_temperature : xr.DataArray
        The ambient temperature data.
    max_forward_temperature : xr.DataArray
        The maximum forward temperature.
    min_forward_temperature : xr.DataArray
        The minimum forward temperature.
    fixed_return_temperature : xr.DataArray
        The fixed return temperature.
    lower_threshold_ambient_temperature : float
        Forward temperature is `max_forward_temperature` for ambient temperatures lower-or-equal this threshold.
    upper_threshold_ambient_temperature : float
        Forward temperature is `min_forward_temperature` for ambient temperatures higher-or-equal this threshold.
    """

    def __init__(
        self,
        ambient_temperature: xr.DataArray,
        max_forward_temperature: float,
        min_forward_temperature: float,
        fixed_return_temperature: float,
        lower_threshold_ambient_temperature: float,
        upper_threshold_ambient_temperature: float,
        rolling_window_ambient_temperature: int,
    ) -> None:
        """
        Initialize the CentralHeatingTemperatureApproximator.

        Parameters
        ----------
        ambient_temperature : xr.DataArray
            The ambient temperature data.
        max_forward_temperature : xr.DataArray
            The maximum forward temperature.
        min_forward_temperature : xr.DataArray
            The minimum forward temperature.
        fixed_return_temperature : xr.DataArray
            The fixed return temperature.
        lower_threshold_ambient_temperature : float
            Forward temperature is `max_forward_temperature` for ambient temperatures lower-or-equal this threshold.
        upper_threshold_ambient_temperature : float
            Forward temperature is `min_forward_temperature` for ambient temperatures higher-or-equal this threshold.
        rolling_window_ambient_temperature : int
            Rolling window size for averaging ambient temperature.
        """

        if any(max_forward_temperature < min_forward_temperature):
            raise ValueError(
                "max_forward_temperature must be greater than min_forward_temperature"
            )
        if any(min_forward_temperature < fixed_return_temperature):
            raise ValueError(
                "min_forward_temperature must be greater than fixed_return_temperature"
            )
        self._ambient_temperature = ambient_temperature
        self.max_forward_temperature = max_forward_temperature
        self.min_forward_temperature = min_forward_temperature
        self.fixed_return_temperature = fixed_return_temperature
        self.lower_threshold_ambient_temperature = lower_threshold_ambient_temperature
        self.upper_threshold_ambient_temperature = upper_threshold_ambient_temperature
        self.rolling_window_ambient_temperature = rolling_window_ambient_temperature

    def ambient_temperature_rolling_mean(self) -> xr.DataArray:
        """
        Property to get ambient temperature.

        Returns
        -------
        xr.DataArray
            Rolling mean of ambient temperature input.
        """
        # bfill to avoid NAs in the beginning
        return (
            self._ambient_temperature.rolling(
                time=self.rolling_window_ambient_temperature
            )
            .mean(skip_na=True)
            .bfill(dim="time")
        )

    @property
    def forward_temperature(self) -> xr.DataArray:
        """
        Property to get dynamic forward temperature.

        Returns
        -------
        xr.DataArray
            Dynamic forward temperatures
        """
        return self._approximate_forward_temperature()

    @property
    def return_temperature(self) -> float:
        """
        Property to get return temperature.

        Returns
        -------
        float
            Return temperature.
        """
        return self._approximate_return_temperature()

    def _approximate_forward_temperature(self) -> xr.DataArray:
        """
        Approximate dynamic forward temperature based on reference curve. Adapted from [Pieper et al. (2019)](https://doi.org/10.1016/j.energy.2019.03.165).

        Returns
        -------
        xr.DataArray
            Dynamic forward temperatures.
        """

        forward_temperature = xr.where(
            self.ambient_temperature_rolling_mean()
            <= self.lower_threshold_ambient_temperature,
            self.max_forward_temperature,
            xr.where(
                self.ambient_temperature_rolling_mean()
                >= self.upper_threshold_ambient_temperature,
                self.min_forward_temperature,
                self.min_forward_temperature
                + (self.max_forward_temperature - self.min_forward_temperature)
                * (
                    self.upper_threshold_ambient_temperature
                    - self.ambient_temperature_rolling_mean()
                )
                / (
                    self.upper_threshold_ambient_temperature
                    - self.lower_threshold_ambient_temperature
                ),
            ),
        )
        return forward_temperature

    def _approximate_return_temperature(self) -> float:
        """
        Approximate return temperature.

        Returns
        -------
        float
            Return temperature.
        """
        return self.fixed_return_temperature

    def _approximate_forward_temperature(self) -> xr.DataArray:
        """
        Approximate dynamic forward temperature.

        Returns
        -------
        xr.DataArray
            Dynamic forward temperatures.
        """
        forward_temperature = xr.where(
            self.ambient_temperature_rolling_mean()
            <= self.lower_threshold_ambient_temperature,
            self.max_forward_temperature,
            xr.where(
                self.ambient_temperature_rolling_mean()
                >= self.upper_threshold_ambient_temperature,
                self.min_forward_temperature,
                self.min_forward_temperature
                + (self.max_forward_temperature - self.min_forward_temperature)
                * (
                    self.upper_threshold_ambient_temperature
                    - self.ambient_temperature_rolling_mean()
                )
                / (
                    self.upper_threshold_ambient_temperature
                    - self.lower_threshold_ambient_temperature
                ),
            ),
        )
        return forward_temperature

    def _approximate_return_temperature(self) -> float:
        """
        Approximate return temperature.

        Returns
        -------
        float
            Return temperature.
        """
        return self.fixed_return_temperature
