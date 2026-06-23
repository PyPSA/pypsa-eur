# supply_curve_v13_ets — EU ETS price sensitivity

## What this is

A one-dimensional sensitivity on the **exogenous EU ETS (ETS-1) carbon price path**
(`costs.emission_prices.co2`), holding everything else identical to **v12 medium**:
same cost CSVs (`custom_costs_medium_hoB.csv`), same gas-for-industry CC growth cap
(1.0 GW/yr), same credit-price sweep (S00–S10, €0–500/t), same solver options.

Only three things differ from the v12 medium configs:
1. `run.name` tag: `t10` → `t10L` (low) / `t10H` (high) / `t10L430` (low, 2050=430) / `t10H525` (high, 2050=525)
2. `run.prefix`: `supply_curve_v12` → `supply_curve_v13_ets`
3. `costs.emission_prices.co2` triple (see below)

Generated reproducibly by `jobs/gen_supply_curve_v13_ets.sh` from the v12 medium configs.

## ETS price paths (EUR/tCO₂, 2030 / 2040 / 2050)

| Path | 2030 | 2040 | 2050 | Stored as | Source / anchor |
|------|-----:|-----:|-----:|-----------|-----------------|
| **low**     |  70 | 130 | 400 | `ets_low/`  (tag `t10L`) | Enerdata POLES core path + removals-sensitivity lower bound |
| **low430**  |  70 | 130 | 430 | `ets_low430/`  (tag `t10L430`) | additional 2050-endpoint point (keeps low 2030/2040) |
| **central** | 119 | 279 | 463 | *(= v12 medium, not regenerated)* | Günther et al. 2025 (PRIMES/PIK), *Climate Policy* — single peer-reviewed trajectory covering all three years |
| **high525** | 160 | 400 | 525 | `ets_high525/` (tag `t10H525`) | additional 2050-endpoint point (keeps high 2030/2040) |
| **high**    | 160 | 400 | 630 | `ets_high/` (tag `t10H`) | LSEG (€160 in 2030; >€400 in 2040 under the 90% target) + Enerdata removals-sensitivity upper bound |

Each year brackets monotonically: 70<119<160, 130<279<400, 400<463<630.
The CAKE (~€460) / PRIMES (€463) agreement in 2050 from two independent models is a
robustness point worth citing. PIK also reports ETS-2 2030 estimates spanning €51–391/t
across studies — the justification for running a *range* rather than a single path.

**low430 (70/130/430) and high525 (160/400/525)** add two extra 2050-endpoint points,
varying only the 2050 value while holding 2030/2040 fixed at the low/high paths. They give
finer resolution of how the binding credit price responds to the 2050 ETS level (the year
where the storage-competition channel bites hardest). Source anchor for the specific 2050
values still to be filled in (currently chosen for extra resolution, not a new citation).

### 2050 infill sweep (pin the binding-credit-price flips)

| Path | 2030 | 2040 | 2050 | Stored as | Purpose |
|------|-----:|-----:|-----:|-----------|---------|
| low405 |  70 | 130 | 405 | `ets_low405/` (tag `t10L405`) | resolve the low-path **100→150 €/t** binding-price flip |
| low420 |  70 | 130 | 420 | `ets_low420/` (tag `t10L420`) | "" |
| low435 |  70 | 130 | 435 | `ets_low435/` (tag `t10L435`) | "" |
| low450 |  70 | 130 | 450 | `ets_low450/` (tag `t10L450`) | "" |
| high500 | 160 | 400 | 500 | `ets_high500/` (tag `t10H500`) | resolve the high-path binding-price flip |
| high550 | 160 | 400 | 550 | `ets_high550/` (tag `t10H550`) | "" |

These hold 2030/2040 fixed and step only the 2050 ETS level, so the binding credit price can
be located to within one infill step rather than relying on the 400→630 endpoints. To find
the flip you do **not** need the full S00–S10 sweep at each point: submit only the credit
prices that bracket the suspected flip, e.g. `S02` (100 €/t) and `S03` (150 €/t) for the low
infill. Whichever is the lowest credit price meeting the policy-goal volume is the binding
price at that ETS level; the flip sits between the two ETS values where that answer changes.

## ⚠️ Interpretation caveats (state these in the chapter)

1. **Standalone CDR is ETS-invariant by construction.** With `cdr_credit_standalone: true`,
   the DAC/biogenic link gets `marginal_cost += co2_price × withdrawal`
   (`scripts/prepare_sector_network.py`), which exactly cancels the ETS revenue from the
   `co2 atmosphere` store. So varying ETS does **not** change credited-CDR profitability.
   What it moves: the fossil generation mix, electricity prices (→ DAC energy opex), and
   **fossil CCS competition for the binding sequestration slots** (fossil CCS avoids ETS,
   so higher ETS makes it bid harder for storage). That contest is the informative result.

2. **ETS here is exogenous; CDR is a separate standalone credit.** External benchmarks where
   CDR *entering* the ETS suppresses the EUA price (ABN AMRO/Bloomberg EUCPM; BNEF) cannot be
   reproduced by this model and belong in the Discussion as comparison points, not inputs.

## How to run

```bash
# (re)generate configs from v12 medium
bash jobs/gen_supply_curve_v13_ets.sh

# submit — respect the ~6-7 concurrent-job HPC limit, so go in waves:
bash jobs/run_supply_curve_v13_ets.sh low       # 11 jobs (S00-S10, low ETS)
bash jobs/run_supply_curve_v13_ets.sh high      # 11 jobs (S00-S10, high ETS)
bash jobs/run_supply_curve_v13_ets.sh low430    # 11 jobs (S00-S10, low ETS 2050=430)
bash jobs/run_supply_curve_v13_ets.sh high525   # 11 jobs (S00-S10, high ETS 2050=525)
# or a subset of credit prices:
bash jobs/run_supply_curve_v13_ets.sh low430 S00 S04 S08

# 2050 infill — only the bracketing credit prices are needed to locate the flip:
bash jobs/run_supply_curve_v13_ets.sh low405 low420 low435 low450 S02 S03   # 8 jobs
bash jobs/run_supply_curve_v13_ets.sh high500 high550 S03 S04               # 4 jobs
```

Results: `results/supply_curve_v13_ets/<run.name>/`. Central comparator is the existing
`results/supply_curve_v12/S*t10hoB-...` runs — no re-run needed.
