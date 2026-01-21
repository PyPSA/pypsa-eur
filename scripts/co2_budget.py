# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def bound_value_for_horizon(bound: Any, current_horizon: int) -> float | None:
    if bound is None:
        return None
    if isinstance(bound, (int, float)):
        return float(bound)
    if isinstance(bound, Mapping):
        if current_horizon in bound:
            value = bound[current_horizon]
            if value is None:
                return None
            return float(value)
        horizon_str = str(current_horizon)
        if horizon_str in bound:
            value = bound[horizon_str]
            if value is None:
                return None
            return float(value)
        return None
    raise TypeError(
        "co2_budget bounds must be null, a number, or a dict mapping year to value. "
        f"Received {type(bound).__name__}."
    )


def co2_budget_for_horizon(
    co2_budget: dict,
    *,
    current_horizon: int,
    baseline_1990: float | None = None,
) -> tuple[float | None, float | None]:
    relative = co2_budget["relative"]

    upper_raw = co2_budget["upper"]
    lower_raw = co2_budget["lower"]

    upper = bound_value_for_horizon(upper_raw, current_horizon)
    lower = bound_value_for_horizon(lower_raw, current_horizon)

    if relative:
        if baseline_1990 is None:
            raise ValueError(
                "co2_budget.relative is true but no 1990 baseline emissions were provided."
            )
        if upper is not None:
            upper *= baseline_1990
        if lower is not None:
            lower *= baseline_1990

    if upper is None:
        # If upper is explicitly disabled but a lower bound is provided, this is invalid.
        if upper_raw is None and lower is not None:
            raise ValueError(
                f"Cannot apply only lower CO2 constraint for horizon {current_horizon}. "
                "The model requires an upper constraint to apply a lower constraint."
            )

        # When no upper bound is configured for this horizon (e.g. missing mapping entry),
        # skip CO2 constraints entirely, regardless of whether a lower bound was provided.
        return None, None

    if lower is not None and upper is not None and lower >= upper:
        raise ValueError(
            f"Lower bound ({lower}) must be less than upper bound ({upper}) "
            f"for horizon {current_horizon}."
        )

    return upper, lower
