# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import xarray as xr

logger = logging.getLogger(__name__)


class PtesTemperatureApproximator:
    """
    A unified class to handle pit thermal energy storage (PTES) temperature-related calculations.

    It calculates top and bottom temperature profiles and approximates storage capacity based on
    temperature differences.

    Attributes
    ----------
    forward_temperature : xr.DataArray
        The forward temperature profile from the district heating network.
    return_temperature : xr.DataArray
        The return temperature profile from the district heating network.
    top_temperature : float | str
        Operational temperature specification for top layer in PTES.
    bottom_temperature : float | str
        Operational temperature specification for bottom layer in PTES.
    temperature_dependent_capacity : bool
        Whether storage capacity varies with temperature. If False, assumes constant capacity.
    design_top_temperature : float
        Maximum design temperature for the top layer of PTES, used for capacity normalization.
    design_bottom_temperature : float
        Minimum design temperature for the bottom layer of PTES, used for capacity normalization.
    """

    def __init__(
        self,
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        top_temperature: float | str,
        bottom_temperature: float | str,
        temperature_dependent_capacity: bool,
        design_top_temperature: float,
        design_bottom_temperature: float,
    ):
        """
        Initialize PtesTemperatureApproximator.

        Parameters
        ----------
        forward_temperature : xr.DataArray
            The forward temperature profile from the district heating network.
        return_temperature : xr.DataArray
            The return temperature profile from the district heating network.
        top_temperature : float | str
            Operational temperature of top layer in PTES. Either a float value or 'forward' for dynamic profiles.
        bottom_temperature : float | str
            Operational temperature of bottom layer in PTES. Either a float value or 'return' for dynamic profiles.
        temperature_dependent_capacity : bool
            Whether storage capacity varies with temperature. If False, assumes constant capacity.
        design_top_temperature : float
            Maximum design temperature for the top layer of PTES, used for capacity normalization
            and clipping dynamic top temperature profiles.
        design_bottom_temperature : float
            Minimum design temperature for the bottom layer of PTES, used for capacity normalization.
        """
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature
        self.top_temperature = top_temperature
        self.bottom_temperature = bottom_temperature
        self.temperature_dependent_capacity = temperature_dependent_capacity
        self.design_top_temperature = design_top_temperature
        self.design_bottom_temperature = design_bottom_temperature

    @property
    def top_temperature_profile(self) -> xr.DataArray:
        """
        PTES top layer temperature profile.

        Returns either the forward temperature (if top_temperature == 'forward'),
        clipped to the design_top_temperature, or a constant temperature profile
        (if top_temperature is a numeric value).

        Returns
        -------
        xr.DataArray
            The resulting top temperature profile for PTES.
        """
        if self.top_temperature == "forward":
            logger.info(
                f"PTES top temperature profile: Using dynamic forward temperature from district heating network, clipped by design top temperature to {self.design_top_temperature}°C "
                f"Forward temperature range: {float(self.forward_temperature.min().values):.1f}°C to {float(self.forward_temperature.max().values):.1f}°C)"
            )
            return self.forward_temperature.where(
                self.forward_temperature <= self.design_top_temperature,
                self.design_top_temperature,
            )
        elif isinstance(self.top_temperature, (int, float)):
            logger.info(
                f"PTES top temperature profile: Using constant temperature of {self.top_temperature}°C "
                f"for all {self.forward_temperature.size} snapshots and nodes"
            )
            return xr.full_like(self.forward_temperature, self.top_temperature)
        else:
            raise ValueError(
                f"Invalid top_temperature: {self.top_temperature}. "
                "Must be 'forward' or a numeric value."
            )

    @property
    def bottom_temperature_profile(self) -> xr.DataArray:
        """
        PTES bottom layer temperature profile.

        Returns either the return temperature (if bottom_temperature == 'return')
        or a constant temperature profile (if bottom_temperature is a numeric value).

        Returns
        -------
        xr.DataArray
            The resulting bottom temperature profile for PTES.
        """
        if self.bottom_temperature == "return":
            logger.info(
                f"PTES bottom temperature profile: Using dynamic return temperature from district heating network "
                f"(shape: {self.return_temperature.shape}, range: {float(self.return_temperature.min().values):.1f}°C to {float(self.return_temperature.max().values):.1f}°C)"
            )
            return self.return_temperature
        elif isinstance(self.bottom_temperature, (int, float)):
            logger.info(
                f"PTES bottom temperature profile: Using constant temperature of {self.bottom_temperature}°C "
                f"for all {self.return_temperature.size} snapshots and nodes"
            )
            return xr.full_like(self.return_temperature, self.bottom_temperature)
        else:
            raise ValueError(
                f"Invalid bottom_temperature: {self.bottom_temperature}. "
                "Must be 'return' or a numeric value."
            )

    @property
    def e_max_pu(self) -> xr.DataArray:
        """
        Calculate e_max_pu for PTES as design_temperature_delta / actual_temperature_delta.

        Returns
        -------
        xr.DataArray
            Normalized delta T values between 0 and 1, representing the
            available storage capacity as a fraction of maximum design capacity.
            If temperature_dependent_capacity is False, returns constant capacity of 1.0.
        """
        if self.temperature_dependent_capacity:
            delta_t = self.top_temperature_profile - self.bottom_temperature_profile
            # Get max possible delta_t for normalization
            max_delta_t = self.design_top_temperature - self.design_bottom_temperature
            normalized_delta_t = delta_t / max_delta_t
            result = normalized_delta_t.clip(min=0)  # Ensure non-negative values
            logger.info(
                f"PTES capacity (e_max_pu): Calculating temperature-dependent capacity. "
                f"Normalization: max_delta_t={max_delta_t:.2f}K (max_top={self.design_top_temperature:.2f}°C, min_bottom={self.design_bottom_temperature:.2f}°C). "
                f"Resulting capacity range: {float(result.min().values):.3f} to {float(result.max().values):.3f} p.u."
            )
            return result
        else:
            logger.info(
                f"PTES capacity (e_max_pu): Using constant capacity of 1.0 p.u. (temperature-independent) "
                f"for all {self.forward_temperature.size} snapshots and nodes"
            )
            return xr.ones_like(self.forward_temperature)
