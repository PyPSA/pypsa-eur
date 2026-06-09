<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Techno-Economic Assumptions

The database of cost assumptions is retrieved from the repository
[PyPSA/technology-data](https://github.com/pypsa/technology-data) and then
saved to a file `data/costs/*/costs_{year}.csv`. The `config/config.yaml` provides options
to choose a reference year. To select a specific version of the cost assumptions, see [Managing Data Versions](data_sources.md#managing_data_versions).

```yaml
{{ yaml_section("costs.year") }}
```

The file includes cost assumptions for all included technologies for specific
years compiled from various sources, namely for

- discount rate,
- lifetime,
- investment (CAPEX),
- fixed operation and maintenance (FOM),
- variable operation and maintenance (VOM),
- fuel costs,
- efficiency, and
- carbon-dioxide intensity.

Many values are taken from a database published by the Danish Energy Agency ([DEA](https://ens.dk/en/our-services/projections-and-models/technology-data)).

The given overnight capital costs are annualised to net present costs
with a discount rate of $r$ over the economic lifetime $n$ using the annuity factor

$$
a = \frac{1-(1+r)^{-n}}{r}.
$$

Based on the parameters above the `marginal_cost` and `capital_cost` of the
system components are automatically calculated.

## Modifying Assumptions

Some cost assumptions (e.g. marginal cost and capital cost) can be directly
set in the `config/config.yaml` (cf. Section [costs](configuration.md#costs_cf) in
[Configuration](configuration.md#config)).
