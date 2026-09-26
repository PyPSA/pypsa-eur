<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Plotting and Summaries

The rules on this page summarise and visualise solved networks and write under
`results/{run}/`.

## Summaries

{{ rules("make_summary", "plot_summary", "plot_base_statistics") }}

## Maps

{{ rules("plot_base_network", "plot_clustered_network", "plot_power_network", "plot_hydrogen_network", "plot_gas_network", "plot_balance_map", "plot_balance_map_interactive", "plot_heat_source_map") }}

## Time series

{{ rules("plot_balance_timeseries", "plot_heatmap_timeseries", "plot_interactive_bus_balance") }}

## Heating diagnostics

{{ rules("build_ambient_air_temperature_yearly_average", "plot_cop_profiles") }}
