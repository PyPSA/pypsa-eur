# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import xarray as xr

logger = logging.getLogger(__name__)


class PtesTemperatureApproximator:
    """
    A unified class to handle pit thermal energy storage (PTES) temperature-related calculations.

    It calculates top temperature profiles, determines when charge or discharge boosting is needed,
    and approximates storage capacity based on temperature differences.

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
    charge_boosting_required : bool
        Whether charge boosting is required/allowed.
    discharge_boosting_required : bool
        Whether discharge boosting is required/allowed.
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
        charge_boosting_required: bool,
        discharge_boosting_required: bool,
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
        charge_boosting_required : bool
            Whether charge boosting is required/allowed.
        discharge_boosting_required : bool
            Whether discharge boosting is required/allowed.
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
        self.charge_boosting_required = charge_boosting_required
        self.discharge_boosting_required = discharge_boosting_required
        self.temperature_dependent_capacity = temperature_dependent_capacity
        self.design_top_temperature = design_top_temperature
        self.design_bottom_temperature = design_bottom_temperature

        if self.charge_boosting_required:
            raise NotImplementedError(
                "Charge boosting for PTES is currently not supported but might be retintroduced in the future."
            )

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

    @property
    def boost_per_discharge(self) -> xr.DataArray:
        """
        Calculate the additional lift required between the store's
        current top temperature and the forward temperature with the lift
        already achieved inside the store.

        Notes
        -----
        The total thermal output required to reach the forward temperature is:

            Q_total = Q_discharge + Q_boost

        The total heat transfer is partitioned into:

            Q_discharge = Ṽ·ρ·cₚ·(T_top − T_bottom)
            Q_boost  = Ṽ·ρ·cₚ·(T_forward − T_top)

        Solving this to constant Ṽ gives α as the ratio of required boost to available store energy:

            α = Q_boost / Q_discharge
              = (T_forward − T_top) / (T_top − T_bottom)

        This expression quantifies the share of PTES output that is covered
        by stored energy relative to the additional heating needed to meet
        the desired forward temperature.

        Returns
        -------
        xr.DataArray
            The resulting fraction of PTES charge that must be further heated.
        """
        if self.discharge_boosting_required:
            result = (
                (self.forward_temperature - self.top_temperature_profile)
                / (self.top_temperature_profile - self.bottom_temperature_profile)
            ).where(self.forward_temperature > self.top_temperature_profile, 0)

            # Count how many snapshots require boosting
            boosting_needed = (
                (self.forward_temperature > self.top_temperature_profile).sum().values
            )
            total_snapshots = self.forward_temperature.size

            logger.info(
                f"Discharge boosting (boost_per_discharge): Enabled. "
                f"Boosting required for {int(boosting_needed)}/{total_snapshots} snapshot-node combinations "
                f"(ratio range: {float(result.min().values):.3f} to {float(result.max().values):.3f})"
            )
            return result
        else:
            logger.info(
                f"Discharge boosting (boost_per_discharge): Not required. "
                f"Returning boost_per_discharge=0 for all {self.forward_temperature.size} snapshots and nodes"
            )
            return xr.zeros_like(self.forward_temperature)

    @property
    def boost_per_charge(self) -> xr.DataArray:
        """
        Calculate the additional boost energy ratio required during charging.

        .. note::
            Charge boosting is currently not implemented and will raise
            NotImplementedError if charge_boosting_required is True.

        This calculates how much additional energy is needed to raise the PTES
        top layer from the forward temperature to the design top temperature,
        relative to the energy delivered by charging.

        Notes
        -----
        To fill the storage from the return temperature all the way up to its
        design top temperature, the total thermal energy required is split into:

            Q_charge = Ṽ·ρ·cₚ·(T_forward − T_bottom)
            Q_boost  = Ṽ·ρ·cₚ·(T_top − T_forward)

        - Q_charge is the energy delivered by charging to the forward setpoint.
        - Q_boost is the extra boost energy needed to reach the design top temperature.

        Defining α as the ratio of boost energy to charge energy:

            α = Q_boost / Q_charge
              = (T_top − T_forward) / (T_forward − T_return)

        Wherever the forward temperature meets or exceeds the design top temperature,
        α is set to zero since no further boost is needed.

        Returns
        -------
        xr.DataArray
            The ratio of additional boost energy needed to the energy delivered by charging.
        """
        if self.charge_boosting_required:
            # Get the max top temperature value
            max_top = (
                self.top_temperature
                if isinstance(self.top_temperature, (int, float))
                else self.top_temperature_profile
            )
            result = (
                (
                    (max_top - self.forward_temperature)
                    / (self.forward_temperature - self.return_temperature)
                )
                .where(self.forward_temperature < max_top, 0)
                .clip(max=1)
            )

            # Count how many snapshots require boosting
            boosting_needed = (self.forward_temperature < max_top).sum().values
            total_snapshots = self.forward_temperature.size

            logger.info(
                f"Charge boosting (boost_per_charge): Enabled. "
                f"Boosting required for {int(boosting_needed)}/{total_snapshots} snapshot-node combinations "
                f"(ratio range: {float(result.min().values):.3f} to {float(result.max().values):.3f})"
            )
            return result
        else:
            logger.info(
                f"Charge boosting (boost_per_charge): Not required. "
                f"Returning boost_per_charge=0 for all {self.forward_temperature.size} snapshots and nodes"
            )
            return xr.zeros_like(self.forward_temperature)
