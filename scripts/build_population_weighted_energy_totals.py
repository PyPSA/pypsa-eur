"""Build population-weighted energy totals."""

import pandas as pd

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_population_weighted_energy_totals',
            simpl='',
            clusters=48,
        )

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    energy_totals = pd.read_csv(snakemake.input.energy_totals, index_col=0)

    nodal_energy_totals = energy_totals.loc[pop_layout.ct].fillna(0.)
    nodal_energy_totals.index = pop_layout.index
    nodal_energy_totals = nodal_energy_totals.multiply(pop_layout.fraction, axis=0)

    nodal_energy_totals.to_csv(snakemake.output[0])
