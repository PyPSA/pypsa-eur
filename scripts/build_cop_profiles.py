"""Build COP time series for air- or ground-sourced heat pumps."""

import xarray as xr


def coefficient_of_performance(delta_T, source='air'):
    """
    COP is function of temp difference source to sink.
    The quadratic regression is based on Staffell et al. (2012)
    https://doi.org/10.1039/C2EE22653G.
    """
    if source == 'air':
        return 6.81 - 0.121 * delta_T + 0.000630 * delta_T**2
    elif source == 'soil':
        return 8.77 - 0.150 * delta_T + 0.000734 * delta_T**2
    else:
        raise NotImplementedError("'source' must be one of  ['air', 'soil']")


if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_cop_profiles',
            simpl='',
            clusters=48,
        )

    for area in ["total", "urban", "rural"]:

        for source in ["air", "soil"]:

            source_T = xr.open_dataarray(
                snakemake.input[f"temp_{source}_{area}"])

            delta_T = snakemake.config['sector']['heat_pump_sink_T'] - source_T

            cop = coefficient_of_performance(delta_T, source)

            cop.to_netcdf(snakemake.output[f"cop_{source}_{area}"])
