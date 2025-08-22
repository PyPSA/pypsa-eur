import logging

import pandas as pd
from _helpers import configure_logging

logger = logging.getLogger(__name__)


def calculate_industrial_production(
    years, scenarios, sector, base_year, base_prod, end_year, end_prod
):
    """
    Generic function to calculate industrial production projections

    Parameters
    ----------
    years : range
        Range of years to project
    scenarios : list
        List of scenarios to model
    sector : str
        Industry sector name
    base_year : int
        Starting year for trend calculation
    base_prod : float
        Production in base year
    end_year : int
        End year for trend calculation
    end_prod : float
        Production in end year
    """
    production = pd.DataFrame(index=years, columns=scenarios)
    trend = abs((end_prod - base_prod) / (end_year - base_year))

    for year in years:
        years_from_end = year - end_year
        for scenario in scenarios:
            if scenario == "regain":
                value = end_prod + (trend * years_from_end)
            elif scenario == "maintain":
                value = end_prod
            elif scenario == "deindustrial":
                value = end_prod - (trend * years_from_end)

            # Ensure no negative values
            production.loc[year, scenario] = max(0, value)

    production["sector"] = sector
    production = production.reset_index()

    # Rename 'index' to 'year' and reorder columns
    production = production.rename(columns={"index": "year"})

    # Reorder columns to put 'sector' after 'year'
    cols = ["year", "sector"] + [
        col for col in production.columns if col not in ["year", "sector"]
    ]
    production = production[cols]

    return production.reset_index()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_policies_projections",
            simpl="",
            clusters="39",
            planning_horizons=2050,
            run="baseline",
        )

    configure_logging(snakemake)

    scenarios = ["regain", "maintain", "deindustrial"]
    years = range(2024, 2051)

    # Production data dictionary
    production_data = {
        "steel": (2000, 186, 2023, 126.219),
        "cement": (2000, 234, 2022, 191),
        "chlorine": (2008, 8.448, 2023, 5.717),
        "ammonia": (2008, 13.832, 2023, 8.959),
        "methanol": (2008, 2.207, 2021, 2.1),
        "hvc": (2008, 60, 2022, 47.2),
        #'hvc': (2008, 17.865 + (12.963 * 28/42) + (2.100 * 28/92),
        #        2023, 10.837 + (9.453 * 28/42) + (1.730 * 28/92))
    }

    # Calculate projections for all sectors
    sector_projections = []
    for sector, (base_year, base_prod, end_year, end_prod) in production_data.items():
        projection = calculate_industrial_production(
            years, scenarios, sector, base_year, base_prod, end_year, end_prod
        )
        sector_projections.append(projection)

    # Combine all projections
    industry_prod_scenarios = pd.concat(sector_projections, ignore_index=True)

    # Save to CSV
    industry_prod_scenarios.to_csv(
        snakemake.output.industry_prod_scenarios, index=False
    )
