# Thesis contributions — code map

This repository is a fork of [PyPSA-Eur](https://github.com/pypsa/pypsa-eur). This
file maps each contribution claimed in the thesis to the exact code that implements
it, so the additions can be distinguished from upstream PyPSA-Eur. It is split into
two parts mirroring the thesis structure: **Methodology** (the modelling machinery
added to the optimisation) and **Results** (the scripts and configs that produce the
thesis figures and scenarios).

**Fork point:** all thesis code is the diff from upstream commit `1fa9a4fb`
(the last PyPSA-Eur commit before the thesis work) to `HEAD`. To reproduce this map:

```bash
git diff 1fa9a4fb HEAD -- scripts/ rules/ Snakefile
```

The canonical scenario configs are `config/Myruns/supply_curve_v12/` (CDR-cost
supply curves: low/medium/high) and `config/Myruns/supply_curve_v13_ets/`
(EU-ETS price sensitivity). Earlier iterations (v6–v11, archives) and intermediate
price sweeps are kept on disk but excluded from the handed-in tree (see `.gitignore`).

---

# Methodology

*Methodological contributions — new modelling machinery added to the PyPSA-Eur
optimisation layer.*

## 1. CDR credit accounting framework
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

**Config interface — [`scripts/lib/validation/config/sector.py`](scripts/lib/validation/config/sector.py):**
`cdr_credit_price`, `cdr_credit_scope`, `cdr_credit_timing`, `cdr_credit_prices_by_scope`.

> **Caveat for the write-up.** Solve-time accounting could fall back to a *capture
> proxy* that may exceed physical geological sequestration. The corrected accounting
> is recomputed post-solve by
> [`regenerate_cdr_accounting.py`](scripts/regenerate_cdr_accounting.py) and
> [`repatch_cdr_accounting.py`](scripts/repatch_cdr_accounting.py). 

---

## 2. Spatially-resolved CO₂ sequestration economics
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

---

## 3. DAC (direct air capture) enhancements
*Country-level weather-dependent performance scaling and configurable technology
variants, beyond upstream's single DAC technology.*

**Code — [`scripts/prepare_sector_network.py`](scripts/prepare_sector_network.py):**
- [`add_dac()`](scripts/prepare_sector_network.py#L1352)
- Weather-factor loading / application: [`_load_dac_weather_factors()`](scripts/prepare_sector_network.py#L1321), [`_node_weather_factors()`](scripts/prepare_sector_network.py#L1335)

**Config interface:** `dac_variants`, `dac_weather_factors` in
[`sector.py`](scripts/lib/validation/config/sector.py).

---

# Results

*Scripts and configs that turn solved networks into the thesis results — supply
curves, decomposition figures, and the EU-ETS sensitivity analysis.*

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

---

## 2. EU-ETS price sensitivity
*Not separate code — an application of the CDR-credit-price machinery
(Methodology §1) driven by configuration.* The ETS price trajectories and binding-price analysis are encoded in
`config/Myruns/supply_curve_v13_ets/` (Low/Central/High ETS paths) and visualised by
the `fig_ets_sensitivity_*` figure scripts.

---

## Scenario configuration (reproducibility)
- Scenario generator: [`config/Myruns/generate_supply_curve_configs.py`](config/Myruns/generate_supply_curve_configs.py)
- Custom optimisation hook: [`data/custom_extra_functionality.py`](data/custom_extra_functionality.py)
- LSF job scripts for the canonical runs: `jobs/run_supply_curve_v12*.sh`, `jobs/run_supply_curve_v13_ets.sh`, `jobs/gen_supply_curve_v13_ets.sh`
</content>
