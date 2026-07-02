# Thesis contributions — code map

This repository is a fork of [PyPSA-Eur](https://github.com/pypsa/pypsa-eur). This
file maps each contribution claimed in the thesis to the exact code that implements
it, so the additions can be distinguished from upstream PyPSA-Eur. It mirrors the
thesis structure:

- **Methodology** — the model extensions added to the optimisation, plus the
  scenario design that drives them.
- **Results** — the scripts that turn solved networks into the thesis figures.

**Fork point:** all thesis code is the diff from upstream commit `1fa9a4fb`
(the last PyPSA-Eur commit before the thesis work) to `HEAD`. To reproduce this map:

```bash
git diff 1fa9a4fb HEAD -- scripts/ rules/ Snakefile
```

The canonical scenario configs live in [`config/Thesis_Runs/`](config/Thesis_Runs):

- [`Scenario_1_Deployment_Response/`](config/Thesis_Runs/Scenario_1_Deployment_Response) —
  CDR deployment response under low / medium / high CDR technology cost.
- [`Robustness_Analysis_ETS/`](config/Thesis_Runs/Robustness_Analysis_ETS) —
  EU-ETS price-sensitivity robustness analysis.

Earlier iterations and intermediate price sweeps are kept on disk (under
`config/Thesis_Runs/01ArchiveSupplyCurve/`) but excluded from the handed-in tree
(see `.gitignore`).

---

# Methodology

## Model extensions

*New modelling machinery added to the PyPSA-Eur optimisation layer.*

### 1. CDR credit accounting framework
*The headline contribution: a new optimisation-layer accounting system for crediting
carbon dioxide removal, with origin-aware eligibility, per-scope price trajectories,
and a choice of crediting timing. Not present in upstream.*

**Core code — [`scripts/solve_network.py`](scripts/solve_network.py) (~+1400 lines):**

| Feature | Location |
|---|---|
| Main credit constraint / objective term | [`add_cdr_credit_accounting()`](scripts/solve_network.py#L647) |
| CO₂ removal origins (dac / biogenic / fossil) | [`CDR_CREDIT_ORIGINS`](scripts/solve_network.py#L415), [`_classify_co2_origin()`](scripts/solve_network.py#L436) |
| Eligibility-aware crediting scope | [`_eligible_cdr_scopes()`](scripts/solve_network.py#L425) |
| Per-scope credit price trajectories | [`_get_cdr_credit_prices_for_period()`](scripts/solve_network.py#L445) |
| Capture vs. sequestration crediting timing | [`add_cdr_credit_accounting()`:651](scripts/solve_network.py#L651) (`cdr_credit_timing`) |
| Physical capture / withdrawal / sequestration tracking | [`_capture_term_data()`](scripts/solve_network.py#L475), [`_withdrawal_term_data()`](scripts/solve_network.py#L539), [`_sequestration_links()`](scripts/solve_network.py#L576) |
| Annual CO₂ sequestration potential limit | [`add_co2_sequestration_limit()`](scripts/solve_network.py#L372) |
| ETS-cancellation surcharge for standalone credit (`cdr_credit_standalone`, applied at network construction) | [`apply_cdr_credit_to_eligible_capture_links()`](scripts/prepare_sector_network.py#L1516) in `prepare_sector_network.py` |

**Config interface — [`scripts/lib/validation/config/sector.py`](scripts/lib/validation/config/sector.py):**
`cdr_credit_price`, `cdr_credit_scope`, `cdr_credit_timing`, `cdr_credit_prices_by_scope`.

**Driven by:** every canonical config — the credit-price sweep `S00`–`S10` sets
`cdr_credit_prices_by_scope` with `cdr_credit_scope: [dac, biogenic]` and
`cdr_credit_timing: sequestration` (see *Scenario design* below).

> **Caveat for the write-up.** Solve-time accounting could fall back to a *capture
> proxy* that may exceed physical geological sequestration. The corrected accounting
> is recomputed post-solve by
> [`regenerate_cdr_accounting.py`](scripts/regenerate_cdr_accounting.py) and
> [`repatch_cdr_accounting.py`](scripts/repatch_cdr_accounting.py).

### 2. Spatially-resolved CO₂ sequestration economics
*Upstream uses a single flat sequestration cost. This contribution adds
country- and node-level cost overrides and location-specific discount rates with
annuity rescaling.*

**Code — [`scripts/prepare_sector_network.py`](scripts/prepare_sector_network.py),
inside [`add_co2_tracking()`](scripts/prepare_sector_network.py#L716):**
- Regional override mechanism: [`apply_regional_overrides()`](scripts/prepare_sector_network.py#L835)
- Sequestration capital cost overrides: [`add_co2_tracking()`:909](scripts/prepare_sector_network.py#L909)
- Node/country discount rates → annuity rescaling: [`add_co2_tracking()`:921](scripts/prepare_sector_network.py#L921)
- CO₂ transport network cost factor: [`add_co2_network()`](scripts/prepare_sector_network.py#L976) (`co2_network_cost_factor`)

**Config interface:** `co2_sequestration_cost_by_country` / `_by_node`,
`co2_sequestration_discount_rate` (+ `_by_country` / `_by_node`),
`co2_network_cost_factor` in [`sector.py`](scripts/lib/validation/config/sector.py).

**Driven by:** the low / medium / high CDR-cost sweep — `co2_sequestration_cost_by_country`
differs across
[`Sensitivity_Low/`](config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_Low),
[`Sensitivity_Medium/`](config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_Medium) and
[`Sensitivity_High/`](config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_High).

### 3. DAC (direct air capture) enhancements
*Country-level weather-dependent performance scaling and configurable technology
variants, beyond upstream's single DAC technology.*

**Code — [`scripts/prepare_sector_network.py`](scripts/prepare_sector_network.py):**
- [`add_dac()`](scripts/prepare_sector_network.py#L1352)
- Weather-factor loading / application: [`_load_dac_weather_factors()`](scripts/prepare_sector_network.py#L1321), [`_node_weather_factors()`](scripts/prepare_sector_network.py#L1335)

**Config interface:** `dac_variants`, `dac_weather_factors` in
[`sector.py`](scripts/lib/validation/config/sector.py).

**Driven by:** every canonical config sets `dac_variants` (DAC-liquidHT / DAC-solidLT)
and `dac_weather_factors` (see *Scenario design* below).

## Scenario design & reproducibility

*The scenarios that exercise the model extensions above, and how to regenerate and
run them. Each scenario family is an eleven-step CDR credit-price sweep (`S00` zero-price
anchor + `S01`–`S10`, currently €0–€500 in €50 steps).*

| Scenario family | Configs |
|---|---|
| CDR deployment response — low CDR cost | [`Sensitivity_Low/`](config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_Low) |
| CDR deployment response — medium CDR cost | [`Sensitivity_Medium/`](config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_Medium) |
| CDR deployment response — high CDR cost | [`Sensitivity_High/`](config/Thesis_Runs/Scenario_1_Deployment_Response/Sensitivity_High) |
| EU-ETS price-sensitivity robustness | [`Robustness_Analysis_ETS/`](config/Thesis_Runs/Robustness_Analysis_ETS) |

**Reproducibility:**
- Scenario generator (regenerates the numbered configs from a template): [`config/Thesis_Runs/generate_supply_curve_configs.py`](config/Thesis_Runs/generate_supply_curve_configs.py)
- Custom optimisation hook: [`data/custom_extra_functionality.py`](data/custom_extra_functionality.py)
- LSF job scripts for the canonical runs: [`jobs/run_scenario1_medium.sh`](jobs/run_scenario1_medium.sh) (+ [`run_scenario1_high.sh`](jobs/run_scenario1_high.sh) / [`run_scenario1_low.sh`](jobs/run_scenario1_low.sh)), [`jobs/run_robustness_ets.sh`](jobs/run_robustness_ets.sh), [`jobs/gen_robustness_ets.sh`](jobs/gen_robustness_ets.sh)

---

# Results

*Scripts that turn the solved networks into the thesis results — supply curves,
decomposition figures, and the EU-ETS sensitivity analysis.*

## 1. Supply-curve construction & CDR post-processing
*Builds the CDR supply curves and decomposition figures that are the thesis results.*

New standalone scripts:
- [`scripts/plot_supply_curve_cdr_split.py`](scripts/plot_supply_curve_cdr_split.py) — supply curve split by CDR origin
- [`scripts/plot_supply_curve_by_technology.py`](scripts/plot_supply_curve_by_technology.py) — supply curve by technology
- [`scripts/plot_cdr_market_gap_storyline.py`](scripts/plot_cdr_market_gap_storyline.py) — CDR market-gap storyline figure

Modified upstream summary scripts (CDR/sequestration accounting in outputs):
[`make_summary.py`](scripts/make_summary.py), [`make_global_summary.py`](scripts/make_global_summary.py).

ETS-sensitivity figures live in `notebooks/fig_ets_sensitivity_*.py` (not tracked in
Git; see project notes).

## 2. EU-ETS price sensitivity
*Not separate code — an application of the CDR-credit-price machinery (Methodology §1)
driven by configuration.* The ETS price trajectories and binding-price analysis are
encoded in the [`Robustness_Analysis_ETS/`](config/Thesis_Runs/Robustness_Analysis_ETS)
configs (Low / Central / High ETS paths; see *Scenario design* above) and visualised by
the `fig_ets_sensitivity_*` figure scripts.
