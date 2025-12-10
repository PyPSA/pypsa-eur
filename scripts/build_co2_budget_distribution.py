# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Build CO2 budget distribution across planning horizons.

This script calculates per-horizon CO2 budget constraints from the unified co2_budget
configuration. It handles:
- Direct per-period caps specified in the configuration
- Total budget distribution using various decay modes (exponential, beta, linear, equal)
- Lower bounds (minimum emissions floors)
- Cumulative total budgets for perfect foresight

The output is a CSV file with columns: horizon, upper_annual, lower_annual, total_cumulative
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import co2_emissions_year

logger = logging.getLogger(__name__)

# CO2 budget configuration constants
VALID_EMISSIONS_SCOPES = frozenset(
    [
        "CO2",
        "CH4",
        "N2O",
        "All greenhouse gases - (CO2 equivalent)",
    ]
)

VALID_DISTRIBUTION_MODES = frozenset(["ex0", "ex1", "be1", "be2", "linear", "equal"])


def _extract_year_values(config_dict: dict) -> dict[int, float]:
    """
    Extract year-based numeric values from config dictionary.

    Parameters
    ----------
    config_dict : dict
        Configuration dictionary that may contain year keys (int or numeric strings)

    Returns
    -------
    dict[int, float]
        Dictionary mapping year (int) to value (float)
    """
    return {
        int(key): float(value)
        for key, value in config_dict.items()
        if isinstance(key, int) or (isinstance(key, str) and key.isdigit())
    }


def parse_co2_budget_config(co2_budget_config: dict) -> dict:
    """
    Parse and normalize CO2 budget configuration to unified structure.

    Parameters
    ----------
    co2_budget_config : dict
        CO2 budget configuration from config file

    Returns
    -------
    dict
        Normalized configuration with keys:
        - emissions_scope: str
        - values: str (interpretation mode: "absolute" or "fraction")
        - upper: dict with keys enable, per_period (dict), total (float), distribution (str)
        - lower: dict with keys enable, per_period (dict)

    Raises
    ------
    KeyError
        If required configuration keys are missing
    """
    # Validate emissions_scope value
    emissions_scope = co2_budget_config["emissions_scope"]
    if emissions_scope not in VALID_EMISSIONS_SCOPES:
        raise ValueError(
            f"Invalid emissions_scope: '{emissions_scope}'. "
            f"Valid options: {', '.join(sorted(VALID_EMISSIONS_SCOPES))}"
        )

    normalized = {
        "emissions_scope": emissions_scope,
        "values": "absolute",  # Default interpretation mode
        "upper": {
            "enable": False,
            "per_period": {},
            "total": None,
            "distribution": None,
        },
        "lower": {"enable": False, "per_period": {}},
    }

    # Get interpretation mode for per-period values (applies to both upper and lower)
    # "absolute": Values are Gt CO2/year (default)
    # "fraction": Values are fractions of 1990 baseline (e.g., 0.450 = 45% of 1990)
    normalized["values"] = co2_budget_config.get("values", "absolute")

    # Parse upper constraints
    if "upper" in co2_budget_config:
        upper_config = co2_budget_config["upper"]

        # FAIL FAST: enable is required when upper section exists
        normalized["upper"]["enable"] = upper_config["enable"]

        # Optional fields for budget distribution (genuinely optional per schema)
        # total: cumulative budget in Gt CO2 (optional - user may only specify per-period caps)
        # distribution: how to split total across periods (optional - only used with total for myopic/overnight)
        # Note: distribution requires total to be set (validated in validate_co2_budget_config)
        normalized["upper"]["total"] = upper_config.get("total")
        normalized["upper"]["distribution"] = upper_config.get("distribution")

        # Extract per-period caps (year keys that are integers)
        per_period_raw = _extract_year_values(upper_config)
        normalized["upper"]["per_period"] = per_period_raw

    # Parse lower constraints
    if "lower" in co2_budget_config:
        lower_config = co2_budget_config["lower"]

        # FAIL FAST: enable is required when lower section exists
        normalized["lower"]["enable"] = lower_config["enable"]

        # Extract per-period floors
        normalized["lower"]["per_period"] = _extract_year_values(lower_config)

    return normalized


def validate_co2_budget_config(config: dict, horizons: list[int]) -> None:
    """
    Validate CO2 budget configuration for consistency and completeness.

    Parameters
    ----------
    config : dict
        Normalized CO2 budget configuration from parse_co2_budget_config
    horizons : list[int]
        Planning horizons

    Raises
    ------
    ValueError
        If configuration is invalid or inconsistent
    """
    upper = config["upper"]
    lower = config["lower"]
    values_mode = config["values"]

    if not upper["enable"]:
        return  # No validation needed if upper constraints disabled

    # Validate values interpretation mode
    if values_mode not in ["absolute", "fraction"]:
        raise ValueError(
            f"Invalid values mode: '{values_mode}'. "
            f"Valid options: 'absolute', 'fraction'"
        )

    # Check that lower < upper where both specified
    # This check works for both fraction and absolute modes since both use the same interpretation
    for year in lower["per_period"]:
        if year in upper["per_period"]:
            if lower["per_period"][year] >= upper["per_period"][year]:
                value_type = "fraction" if values_mode == "fraction" else "Gt CO2/year"
                raise ValueError(
                    f"Lower bound ({lower['per_period'][year]} {value_type}) must be less than "
                    f"upper bound ({upper['per_period'][year]} {value_type}) for year {year}"
                )

    # Validate distribution requirements
    has_total = upper["total"] is not None
    has_distribution = upper["distribution"] is not None
    has_all_period_caps = all(h in upper["per_period"] for h in horizons)

    if has_total and not has_all_period_caps:
        # Total specified but not all period caps provided
        # Distribution is optional (null means cumulative for perfect foresight)
        if has_distribution and upper["distribution"] not in VALID_DISTRIBUTION_MODES:
            raise ValueError(
                f"Invalid distribution mode: '{upper['distribution']}'. "
                f"Valid options: {', '.join(sorted(VALID_DISTRIBUTION_MODES))}"
            )

    if has_distribution and not has_total:
        logger.warning(
            f"Distribution mode '{upper['distribution']}' specified but no total budget provided. "
            "Distribution will be ignored."
        )


def calculate_budget_distribution(
    total_budget: float,
    distribution_mode: str,
    horizons: list[int],
    e_0: float,
) -> dict[int, float]:
    """
    Calculate distributed CO2 budget fractions across planning horizons.

    Parameters
    ----------
    total_budget : float
        Total CO2 budget in Gt CO2
    distribution_mode : str
        Distribution mode (ex0, ex1, be1, be2, linear, equal)
    horizons : list[int]
        Planning horizon years
    e_0 : float
        Initial emissions (Gt CO2), typically from most recent year (e.g., 2018)

    Returns
    -------
    dict[int, float]
        Mapping of horizon year to distributed budget (Gt CO2/year)

    Raises
    ------
    ValueError
        If inputs are invalid or distribution mode is unsupported
    """
    from scipy.stats import beta

    # FAIL FAST: Validate inputs
    if e_0 <= 0:
        raise ValueError(
            f"Initial emissions e_0 must be positive, got {e_0}. "
            "Cannot distribute budget with zero or negative baseline emissions."
        )

    if total_budget < 0:
        raise ValueError(f"Total budget must be non-negative, got {total_budget}")

    if len(horizons) < 1:
        raise ValueError(
            "At least one planning horizon required for budget distribution"
        )

    if distribution_mode == "linear" and len(horizons) == 1:
        raise ValueError(
            "Linear distribution requires at least 2 planning horizons. "
            f"Got {len(horizons)} horizon: {horizons}"
        )

    t_0 = horizons[0]

    if distribution_mode.startswith("ex"):
        # Exponential decay
        r = float(distribution_mode[2:]) if len(distribution_mode) > 2 else 0.0
        T = total_budget / e_0
        m = (1 + np.sqrt(1 + r * T)) / T

        def exponential_decay(t):
            return e_0 * (1 + (m + r) * (t - t_0)) * np.exp(-m * (t - t_0))

        distribution = {t: exponential_decay(t) for t in horizons}

    elif distribution_mode.startswith("be"):
        # Beta decay
        be_param = float(distribution_mode[2:]) if len(distribution_mode) > 2 else 1.0
        t_f = t_0 + (2 * total_budget / e_0).round(0)

        def beta_decay(t):
            cdf_term = (t - t_0) / (t_f - t_0)
            return e_0 * (1 - beta.cdf(cdf_term, be_param, be_param))

        distribution = {t: beta_decay(t) for t in horizons}

    elif distribution_mode == "linear":
        # Linear decrease from e_0 to 0
        n_periods = len(horizons)
        distribution = {
            horizons[i]: e_0 * (1 - i / (n_periods - 1)) for i in range(n_periods)
        }

    elif distribution_mode == "equal":
        # Equal distribution across periods
        per_period = total_budget / len(horizons)
        distribution = {t: per_period for t in horizons}

    else:
        raise ValueError(
            f"Unknown distribution mode: {distribution_mode}. "
            f"Valid options: {', '.join(sorted(VALID_DISTRIBUTION_MODES))}"
        )

    return distribution


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_co2_budget_distribution",
            run="test-run",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Get configuration
    co2_budget_config = snakemake.params.co2_budget
    countries = snakemake.params.countries
    planning_horizons = snakemake.params.planning_horizons
    sector_opts = snakemake.params.sector

    # Parse and validate configuration
    config_normalized = parse_co2_budget_config(co2_budget_config)

    # Handle both single value and list for planning_horizons
    if isinstance(planning_horizons, (int, str)):
        horizons = [int(planning_horizons)]
    else:
        horizons = [int(h) for h in planning_horizons]

    validate_co2_budget_config(config_normalized, horizons)

    upper = config_normalized["upper"]
    lower = config_normalized["lower"]
    values_mode = config_normalized["values"]

    # Initialize output data structure
    distribution_data = []

    # ALWAYS calculate baseline emissions (needed for fractions AND distribution)
    emissions_scope = config_normalized["emissions_scope"]
    e_1990 = None
    e_0 = None
    distributed_budgets = {}

    # Calculate baselines if any constraint is enabled
    if upper["enable"] or lower["enable"]:
        logger.debug(f"Calculating baseline emissions (scope: {emissions_scope})")

        e_1990 = co2_emissions_year(
            countries,
            snakemake.input.eurostat,
            sector_opts,
            emissions_scope,
            snakemake.input.co2,
            year=1990,
        )

        e_0 = co2_emissions_year(
            countries,
            snakemake.input.eurostat,
            sector_opts,
            emissions_scope,
            snakemake.input.co2,
            year=2018,
        )

        logger.debug(
            f"Baseline emissions: 1990={e_1990:.3f} Gt CO2, 2018={e_0:.3f} Gt CO2"
        )

    # Calculate distributed budgets if distribution mode is specified
    if upper["enable"] and upper["total"] and upper["distribution"]:
        distributed_budgets = calculate_budget_distribution(
            upper["total"],
            upper["distribution"],
            horizons,
            e_0,
        )

        logger.debug(f"Distributed budget using mode '{upper['distribution']}':")
        for h, budget in distributed_budgets.items():
            logger.debug(f"  {h}: {budget:.3f} Gt CO2/year")

    # Build distribution data for each horizon
    for horizon in horizons:
        row = {
            "horizon": horizon,
            "upper_annual": np.nan,
            "lower_annual": np.nan,
            "total_cumulative": np.nan,
        }

        if upper["enable"]:
            # Priority 1: Direct per-period cap
            if horizon in upper["per_period"]:
                direct_cap_raw = upper["per_period"][horizon]

                # Convert fraction to absolute value if needed
                if values_mode == "fraction":
                    # Value is a fraction of 1990 emissions (e.g., 0.450 = 45% of 1990)
                    direct_cap = direct_cap_raw * e_1990
                    logger.debug(
                        f"Horizon {horizon}: fraction {direct_cap_raw:.3f} × 1990 baseline "
                        f"{e_1990:.3f} Gt = {direct_cap:.3f} Gt CO2/year"
                    )
                else:
                    # Value is already absolute (Gt CO2/year)
                    direct_cap = direct_cap_raw
                    logger.debug(
                        f"Horizon {horizon}: absolute cap {direct_cap:.3f} Gt CO2/year"
                    )

                # Priority 2: Distributed budget (if available)
                if horizon in distributed_budgets:
                    distributed_cap = distributed_budgets[horizon]
                    # Take minimum of direct and distributed
                    row["upper_annual"] = min(direct_cap, distributed_cap)
                    logger.debug(
                        f"Horizon {horizon}: using min of direct cap ({direct_cap:.3f}) "
                        f"and distributed ({distributed_cap:.3f}) = {row['upper_annual']:.3f} Gt CO2/year"
                    )
                else:
                    row["upper_annual"] = direct_cap

            # No direct cap, use distributed budget if available
            elif horizon in distributed_budgets:
                row["upper_annual"] = distributed_budgets[horizon]
                logger.debug(
                    f"Horizon {horizon}: using distributed budget {row['upper_annual']:.3f} Gt CO2/year"
                )

            # Set total_cumulative only on last horizon if distribution is null
            if horizon == horizons[-1] and upper["total"] and not upper["distribution"]:
                row["total_cumulative"] = upper["total"]
                logger.debug(
                    f"Horizon {horizon}: total cumulative budget {upper['total']:.3f} Gt CO2 "
                    "(for perfect foresight)"
                )

        if lower["enable"] and horizon in lower["per_period"]:
            lower_raw = lower["per_period"][horizon]

            # Convert fraction to absolute value if needed (same logic as upper bounds)
            # Note: We use the top-level "values" setting for both upper and lower bounds
            # to ensure consistent interpretation (mixing fractions and absolutes would be confusing)
            if values_mode == "fraction":
                # Value is a fraction of 1990 emissions
                row["lower_annual"] = lower_raw * e_1990
                logger.debug(
                    f"Horizon {horizon}: lower fraction {lower_raw:.3f} × 1990 baseline "
                    f"{e_1990:.3f} Gt = {row['lower_annual']:.3f} Gt CO2/year"
                )
            else:
                # Value is already absolute (Gt CO2/year)
                row["lower_annual"] = lower_raw
                logger.debug(
                    f"Horizon {horizon}: lower bound {row['lower_annual']:.3f} Gt CO2/year"
                )

        distribution_data.append(row)

    # Create DataFrame with explicit index
    df = pd.DataFrame(distribution_data)
    df = df.set_index("horizon")

    # Validate that at least some constraints were specified
    if upper["enable"]:
        specified_upper = df["upper_annual"].notna().sum()
        if specified_upper == 0:
            logger.warning(
                "co2_budget.upper is enabled but no upper constraints were specified "
                "for any planning horizon. No upper CO2 constraints will be applied."
            )
        else:
            logger.debug(
                f"Upper CO2 constraints specified for {specified_upper}/{len(horizons)} horizons"
            )

    if lower["enable"]:
        specified_lower = df["lower_annual"].notna().sum()
        if specified_lower == 0:
            logger.warning(
                "co2_budget.lower is enabled but no lower constraints were specified "
                "for any planning horizon. No lower CO2 constraints will be applied."
            )
        else:
            logger.debug(
                f"Lower CO2 constraints specified for {specified_lower}/{len(horizons)} horizons"
            )

    # Save to CSV
    df.to_csv(snakemake.output.co2_budget_distribution)
