# supply_curve_v14_netzero — economy-wide net-zero sensitivity

## What this is

A robustness check that turns the **economy-wide CO2 cap back on**. The main analysis
runs with `co2_budget: null` (no cap), so engineered CDR deploys purely on the standalone
credit economics. This sensitivity asks whether that picture survives when the system is
*also* required to decarbonise to net-zero, holding everything else identical to **v12
medium** (cost CSVs `custom_costs_medium_hoB.csv`, central ETS 119/279/463, gas-for-industry
CC cap 1.0, credit sweep S00-S10, solver options).

Only three things differ from the v12 medium configs:
1. `run.name` tag: `t10hoB` -> `t10NZhoB`
2. `run.prefix`: `supply_curve_v12` -> `supply_curve_v14_netzero`
3. `co2_budget`: `null` -> `{2030: 0.45, 2040: 0.10, 2050: 0.0}`

Generated reproducibly by `jobs/gen_robustness_netzero.sh` from the v12 medium configs.

## CO2 budget trajectory (fraction of 1990 emissions)

| Year | Cap (frac. 1990) |
|-----:|-----------------:|
| 2030 | 0.45 |
| 2040 | 0.10 |
| 2050 | 0.00 (net-zero) |

PyPSA-Eur default decarbonisation path. Implemented through `add_co2limit()` as a global
`co2_atmosphere` `GlobalConstraint` per planning year (`scripts/prepare_sector_network.py`).
Myopic foresight: each horizon solves under its own cap, inheriting the previous solve as
brownfield.

## Interpretation caveats (state these in the chapter)

1. **Two carbon signals now coexist.** With a binding cap the model derives an *endogenous*
   CO2 shadow price on top of the exogenous ETS path. The standalone-credit cancellation
   term (`marginal_cost += co2_price x withdrawal`) still uses the **exogenous** ETS price,
   so the credit cancels the exogenous ETS value, not the shadow price. Read the credit
   price as the standalone signal *on top of* whatever the cap enforces.

2. **The cap can pull CDR in by itself.** Under net-zero, DAC/BECCS may deploy to relax the
   atmosphere constraint even at low credit price, so the binding credit price under the cap
   is not directly comparable to the no-cap binding price. The informative result is how the
   supply curve and the fossil/process-CCS competition for storage shift once net-zero is
   imposed.

3. **Feasibility watch.** The 2050 net-zero cap (0.0) plus the 600 Mt storage potential and
   ENSPRESO biomass ceiling could bind hard or, in 2040 (0.10), risk infeasibility. Check
   each solve log; if 2040 is infeasible, the looser `0.45 / 0.25 / 0.0` trajectory is the
   fallback.

## How to run

```bash
# (re)generate configs from v12 medium
bash jobs/gen_robustness_netzero.sh

# submit (respect ~6-7 concurrent-job HPC limit):
bash jobs/run_robustness_netzero.sh                 # all 11 (S00-S10)
bash jobs/run_robustness_netzero.sh S00 S02 S03 S04 # bracketing prices only
```

Results: `results/supply_curve_v14_netzero/<run.name>/`. Comparator is the existing
`results/supply_curve_v12/S*t10hoB-...` runs (same configs, `co2_budget: null`) — no re-run
of the comparator needed.
