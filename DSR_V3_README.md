# DSR v3: Optional DSR with Economic Costs

## Overview

This branch (`industry-dsr-v3-optional`) implements **economically optional DSR** that allows the optimizer to decide whether building DSR capacity is cost-effective.

**Worktree location**: `~/.cursor/worktrees/pypsa-eurelectric__WSL__Ubuntu_/dsr-v3`

---

## Key Changes from v2

### 1. **DSR Capacity is Now Extendable** ✅

**Before (v2)**:
```python
e_nom=e_nom,  # Fixed capacity - always built
# No e_nom_extendable
# No capital_cost
```

**After (v3)**:
```python
e_nom=0.0,  # Start with 0
e_nom_extendable=True,  # Let optimizer decide!
e_nom_max=e_nom,  # Use calculated capacity as upper limit
capital_cost=store_capital_cost,  # Add economic cost
```

**Impact**: Optimizer can now choose whether to build DSR based on cost-benefit analysis.

---

### 2. **Added Economic Parameters** 💰

New config parameters in `config/config.yaml`:

```yaml
industry:
  dsr:
    # DSR ECONOMIC PARAMETERS (NEW in v3)
    store_capital_cost: 30.0  # EUR/MWh/year (energy capacity)
    link_capital_cost: 0.0    # EUR/MW/year (power capacity)
    charge_efficiency: 0.98   # 2% loss per direction
    discharge_efficiency: 0.98
    # Net roundtrip: 96.04%
```

**Rationale**:
- **30 EUR/MWh/year**: Conservative estimate for industrial DSR
  - Control systems, sensors, scheduling software
  - Process modifications if needed
  - Amortized over lifetime
- **2% efficiency loss**: Represents opportunity costs of production scheduling
- **Link cost = 0**: Links are logical connections, minimal physical infrastructure

---

### 3. **Relaxed Checkpoint Frequency** ⏰

**Before (v2)**:
- Checkpoints at **6am AND 6pm** (twice daily)
- Only **12-hour shifting windows**
- 728 checkpoints per year

**After (v3)**:
- Checkpoints at **midnight only** (once daily)
- **24-hour shifting windows**
- 365 checkpoints per year

**Example**:
```yaml
"Iron & steel industry|Scrap-EAF": [0]  # Daily checkpoint at midnight
# Instead of [6, 18]
```

**Benefits**:
- Better captures day/night price arbitrage
- More flexibility for load shifting
- Still ensures daily production completion
- Fewer forced cycling events

---

### 4. **Added Efficiency Losses** 📉

Links now have realistic efficiency parameters:

```python
n.add("Link", charge_links,
      efficiency=0.98)  # 2% loss when charging

n.add("Link", discharge_links,
      efficiency=0.98)  # 2% loss when discharging
```

**Roundtrip efficiency**: 0.98 × 0.98 = **96.04%**

This models the **opportunity costs** of production scheduling flexibility.

---

## How It Works

### Economic Optimization

The model now solves:

```
Minimize: system_cost + DSR_capital_cost + efficiency_losses

Subject to:
  - All normal constraints (generation, transmission, etc.)
  - DSR checkpoints (must return to 0 daily at midnight)
  - DSR capacity limits: 0 ≤ e_nom_opt ≤ e_nom_max
```

**Decision**:
- If `arbitrage_value > (capital_cost + cycling_costs)` → Build DSR
- If `arbitrage_value < costs` → Build zero or minimal DSR
- Optimizer chooses optimal capacity between 0 and maximum

---

### Checkpoint Philosophy

**V3 maintains production discipline** while adding economic optionality:

1. **When DSR is built**: Checkpoints ensure daily production completion
2. **Economic optionality**: DSR is only built if cost-effective
3. **Operational flexibility**: 24-hour windows allow better arbitrage

**Key insight**: Hard checkpoints are OK because they represent **physical production constraints**. The issue wasn't checkpoints themselves, but that v2 **forced building** DSR capacity regardless of economics.

---

## Expected Results

### Scenario A: DSR is Economically Viable

If price volatility justifies DSR costs:

- ✅ **Lower system costs** (DSR provides arbitrage value > its costs)
- ✅ **Less peaking capacity** (OCGT, CCGT)
- ✅ **Lower average prices** (reduced capacity investment)
- ✅ **Higher DSR utilization** (>10-20% average fill level)
- ✅ **Reasonable DSR capacity** (optimizer builds optimal amount)

### Scenario B: DSR is NOT Economically Viable

If costs exceed benefits:

- ✅ **Zero or minimal DSR capacity** (optimizer chooses not to build)
- ✅ **Results similar to baseline** (no DSR impact)
- ✅ **This is OK!** Confirms DSR isn't cost-effective for these parameters

### Comparison to V2

| Metric | V2 (Mandatory) | V3 (Optional) |
|--------|----------------|---------------|
| **DSR Capacity** | 63.62 GWh (fixed) | 0 to 63.62 GWh (optimized) |
| **System Cost Change** | +3.02 B EUR ❌ | TBD (likely reduced or neutral) |
| **Price Change** | +3.74 EUR/MWh ❌ | TBD (likely reduced or neutral) |
| **DSR Utilization** | 1.9% (forced cycling) | TBD (efficient if built) |
| **Peaking Capacity** | +50 GW OCGT ❌ | TBD (likely less increase) |
| **Checkpoints** | 728/year (12h windows) | 365/year (24h windows) |

---

## Configuration File

**Location**: `~/.cursor/worktrees/pypsa-eurelectric__WSL__Ubuntu_/dsr-v3/config/config.yaml`

**Key settings**:
```yaml
run:
  name: "EU_test_run_dsr_v3_optional"

industry:
  dsr:
    enable: true
    store_capital_cost: 30.0  # Can adjust for sensitivity analysis
    charge_efficiency: 0.98
    discharge_efficiency: 0.98
    
    restriction_time:
      # All technologies: once-daily checkpoint at midnight
      "Iron & steel industry|Scrap-EAF": [0]
      # ... etc
```

---

## Running the Model

### Option 1: Quick Test (If Baseline Data Exists)

```bash
cd ~/.cursor/worktrees/pypsa-eurelectric__WSL__Ubuntu_/dsr-v3

# Run solve network only (if base network exists)
snakemake -call results/EU_test_run_dsr_v3_optional/networks/base_s_38___2030.nc \
  --cores all --configfile config/config.yaml
```

### Option 2: Full Run

```bash
cd ~/.cursor/worktrees/pypsa-eurelectric__WSL__Ubuntu_/dsr-v3

# Full workflow
snakemake -call results/EU_test_run_dsr_v3_optional/networks/elec_s_38___2030.nc \
  --cores all --configfile config/config.yaml
```

---

## Parameter Sensitivity Analysis

You can test different cost assumptions by modifying `config.yaml`:

### Conservative (High Cost) - Tests if DSR survives high costs:
```yaml
store_capital_cost: 50.0  # EUR/MWh/year
charge_efficiency: 0.95
discharge_efficiency: 0.95
```

### Aggressive (Low Cost) - Tests maximum DSR potential:
```yaml
store_capital_cost: 10.0  # EUR/MWh/year
charge_efficiency: 0.99
discharge_efficiency: 0.99
```

### Weekly Checkpoints - Tests longer shifting windows:
```yaml
restriction_time:
  "Iron & steel industry|Scrap-EAF": [0]  # Monday at midnight
  # And change snapshots to weekly resolution (168h windows)
```

---

## Comparison Strategy

To validate v3 improvements, compare three runs:

1. **Baseline** (no DSR): `~/work/pypsa-eurelectric/results/EU_test_run`
2. **V2** (mandatory DSR): `~/.cursor/worktrees/kse/results/EU_test_run_dsr`
3. **V3** (optional DSR): `~/.cursor/worktrees/dsr-v3/results/EU_test_run_dsr_v3_optional`

**Key metrics**:
- System cost (objective value)
- DSR capacity built (e_nom_opt)
- DSR utilization (average fill %)
- Average electricity prices
- Peaking capacity (OCGT, CCGT)
- Total curtailment

**Expected outcome**:
- If v3 builds DSR: Cost should be lower than baseline (DSR is beneficial)
- If v3 doesn't build DSR: Cost should match baseline (DSR isn't beneficial)
- Either way: v3 should NOT be worse than baseline (unlike v2)

---

## Code Changes Summary

**Modified files**:
1. `config/config.yaml`:
   - Changed run name to `EU_test_run_dsr_v3_optional`
   - Added `store_capital_cost`, `link_capital_cost`, efficiencies
   - Relaxed checkpoints from twice-daily to once-daily

2. `scripts/prepare_sector_network.py` (~line 5397-5448):
   - Read economic parameters from config
   - Set `e_nom=0.0` and `e_nom_extendable=True`
   - Set `e_nom_max=e_nom` (calculated capacity as upper limit)
   - Added `capital_cost` to stores
   - Added `efficiency` and `capital_cost` to links

**No changes needed**:
- `scripts/solve_network.py` (redundant constraint already removed in v2)
- `scripts/build_industry_dsr_profile.py` (checkpoint generation works as-is)

---

## Debugging Tips

If v3 still shows issues:

1. **Check DSR capacity built**:
   ```python
   dsr_stores = n.stores[n.stores.carrier == 'industry dsr']
   print(f"Total capacity built: {dsr_stores.e_nom_opt.sum():.2f} MWh")
   print(f"Max allowed: {dsr_stores.e_nom_max.sum():.2f} MWh")
   ```

2. **Check capital costs are applied**:
   ```python
   print(f"Capital cost: {dsr_stores.capital_cost.mean():.2f} EUR/MWh/year")
   ```

3. **Check utilization**:
   ```python
   fill_level = n.stores_t.e[dsr_stores.index].mean() / dsr_stores.e_nom_opt
   print(f"Average fill: {fill_level.mean()*100:.1f}%")
   ```

4. **Compare to baseline**:
   - If DSR capacity ≈ 0 → DSR isn't economic, which is OK
   - If DSR capacity > 0 but costs still higher → Check if checkpoint constraints too tight

---

## Next Steps

1. ✅ Config updated with v3 parameters
2. ✅ Code updated in prepare_sector_network.py
3. ⏳ **Run the model** and check results
4. ⏳ **Compare** v3 vs v2 vs baseline
5. ⏳ If needed: **Sensitivity analysis** with different cost assumptions

---

## Branch Management

- **V2 branch** (current working): `industry-dsr-v2` in `~/.cursor/worktrees/kse`
  - Keep as checkpoint - this is the "working but suboptimal" version
  
- **V3 branch** (improvements): `industry-dsr-v3-optional` in `~/.cursor/worktrees/dsr-v3`
  - Test ground for optional DSR improvements
  
- **Baseline**: `master` branch in `~/work/pypsa-eurelectric`
  - Reference for comparison

**Safe to experiment!** V2 is preserved if v3 doesn't work out.
