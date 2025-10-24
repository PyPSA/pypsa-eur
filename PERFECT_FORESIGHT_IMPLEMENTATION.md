<!--
SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

SPDX-License-Identifier: CC-BY-4.0
-->

# Perfect Foresight Implementation in Streamlined Workflow

This document describes the three critical features implemented in `scripts/compose_network.py` to enable perfect foresight optimization in the streamlined PyPSA-EUR workflow.

## Overview

These features were previously only available in the deprecated `scripts/prepare_perfect_foresight.py` and are now integrated into the streamlined workflow's `compose_network.py` script.

## Implemented Features

### 1. Investment Period Objective Weightings with Social Discounting

**Location**: Lines 105-156, 499-549, 1700-1711 in `scripts/compose_network.py`

**Purpose**: Properly weight costs across different investment periods using social discounting to enable fair comparison of investments made at different times.

**Implementation**:

```python
def get_investment_weighting(time_weighting: pd.Series, r: float = 0.01) -> pd.Series:
    """Calculate investment period objective weightings based on social discounting."""
```

**Key Functions**:
- `get_investment_weighting()`: Calculates social discount factors for each period
- `apply_investment_period_weightings()`: Applies both time weightings (years) and objective weightings to the network

**When Applied**: After network concatenation in perfect foresight mode, only when `n._multi_invest` is True

**Configuration**: Uses `params.costs["social_discountrate"]` from config (default: 0.02 or 2%)

**Impact**: Without this, investment decisions in different periods would not be properly discounted, leading to incorrect optimization results.

**Testing**:
```bash
# Verify the function works
python -c "
import pandas as pd
from scripts.compose_network import get_investment_weighting

periods = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
weights = get_investment_weighting(periods, r=0.02)
print('Time periods:', periods.to_dict())
print('Objective weights:', weights.to_dict())
"
```

---

### 2. Store Cycling Adjustments for Multi-Period Optimization

**Location**: Lines 552-620, 1709-1711 in `scripts/compose_network.py`

**Purpose**: Adjust storage cycling behavior for perfect foresight to ensure stores behave correctly across investment periods.

**Implementation**:

```python
def adjust_stores_for_perfect_foresight(n: pypsa.Network) -> None:
    """Adjust store cycling behavior for perfect foresight multi-period optimization."""
```

**Store Behavior Changes**:

1. **Cyclic stores** (e.g., batteries, hydrogen storage):
   - Changed from `e_cyclic=True` to `e_cyclic_per_period=True`
   - Ensures they cycle within each investment period, not over entire horizon

2. **Non-cyclic stores** (CO2, biomass, biogas, EV battery):
   - Set `e_cyclic=False` and `e_cyclic_per_period=False`
   - These accumulate/deplete over the entire planning horizon

3. **Biomass stores**:
   - Set `e_initial_per_period=True`
   - Annual biomass availability is renewed at the start of each period

**When Applied**: After network concatenation in perfect foresight mode, only when `n._multi_invest` is True

**Impact**: Without this, storage constraints would be incorrect - stores would either cycle over the entire 20-30 year horizon (wrong) or not reset annually for renewable resources like biomass (also wrong).

**Affected Stores**:
- Cyclic: H2 storage, batteries, thermal storage, etc.
- Non-cyclic CO2: CO2 capture/storage stores
- Non-cyclic renewable: solid biomass, biogas
- EV batteries: daily patterns, not multi-year cyclic

---

### 3. Multi-Period Time Segmentation

**Location**: Lines 1282-1386, 1479-1505 in `scripts/compose_network.py`

**Purpose**: Apply tsam time segmentation separately for each investment period to preserve multi-period structure.

**Implementation**:

```python
def apply_time_segmentation_multiperiod(
    n: pypsa.Network, segments: int, solver_name: str = "cbc"
) -> None:
    """Apply time segmentation to multi-period networks."""
```

**Key Differences from Single-Period**:
- Processes each investment period separately (preserves period structure)
- Handles MultiIndex snapshots (period, timestep)
- Maintains investment period boundaries

**Auto-Detection**: The code automatically detects if the network has MultiIndex snapshots and uses the appropriate segmentation method:

```python
if isinstance(n.snapshots, pd.MultiIndex):
    # Use multi-period segmentation
    apply_time_segmentation_multiperiod(n, segments, solver_name)
else:
    # Use single-period segmentation
    apply_time_segmentation(n, segments)
```

**When Applied**: During network preparation, if time segmentation is enabled in config

**Impact**: Without this, time segmentation would fail on multi-period networks or collapse all periods into a single timeline.

---

## Integration Flow

For perfect foresight optimization, the features are applied in this order:

```
1. Compose network for each horizon individually
2. Concatenate horizons together (if not first horizon)
   ↓
3. Apply investment period weightings (TASK 1)
   ↓
4. Adjust store cycling behavior (TASK 2)
   ↓
5. Apply multi-period time segmentation if enabled (TASK 3)
   ↓
6. Export network for solving
```

## Configuration

All features use existing configuration parameters:

```yaml
# config/config.yaml
costs:
  social_discountrate: 0.02  # 2% social discount rate

clustering:
  temporal:
    time_segmentation:
      enable: true
      segments: 10  # Number of segments per period
```

## Testing

### Unit Test (Manual)

```bash
# Test compilation
python -m compileall scripts/compose_network.py

# Test individual functions
python -c "
import pandas as pd
import pypsa
from scripts.compose_network import get_investment_weighting, adjust_stores_for_perfect_foresight

# Test investment weighting
periods = pd.Series([10, 10, 10], index=[2030, 2040, 2050])
weights = get_investment_weighting(periods, r=0.02)
assert len(weights) == 3
assert weights[2030] > weights[2050]  # Earlier periods should have higher weight
print('Investment weighting test: PASS')

# Test store adjustment
n = pypsa.Network()
n.add('Bus', 'bus')
n.add('Store', 'battery', bus='bus', carrier='battery', e_cyclic=True)
n.add('Store', 'co2', bus='bus', carrier='co2', e_cyclic=False)
adjust_stores_for_perfect_foresight(n)
assert n.stores.loc['battery', 'e_cyclic_per_period'] == True
assert n.stores.loc['co2', 'e_cyclic_per_period'] == False
print('Store adjustment test: PASS')
"
```

### Integration Test

```bash
# Test electricity-only perfect foresight (when implemented)
snakemake solve_networks --configfile config/test/config.electricity-perfect.yaml --cores 4

# Test sector-coupled perfect foresight (when implemented)
snakemake solve_networks --configfile config/test/config.perfect.yaml --cores 4
```

## Differences from prepare_perfect_foresight.py

### Implemented in This PR

1. ✅ Investment period objective weightings
2. ✅ Store cycling adjustments
3. ✅ Multi-period time segmentation

### NOT Yet Implemented

These features from `prepare_perfect_foresight.py` are NOT yet implemented and may need to be added later:

1. **Phase-out constraints** (`set_all_phase_outs()`):
   - Enforces national energy policy phase-outs (nuclear, coal, lignite)
   - IMPACT: Policy constraints won't be enforced
   - DECISION: Could be added as separate enhancement

2. **HVDC transport model** (`hvdc_transport_model()`):
   - Converts AC lines to DC links for multi-decade optimization
   - IMPACT: Line expansion in perfect foresight may not work correctly
   - DECISION: Needs investigation - may be required for perfect foresight

3. **Heat pump efficiency updates for all periods** (`update_heat_pump_efficiency()`):
   - Different from brownfield version - updates ALL heat pumps to each period's technology
   - IMPACT: Heat pump efficiency may not be handled correctly across periods
   - DECISION: Current brownfield implementation may need modification

4. **H2 boiler retrofits** (`add_H2_boilers()`):
   - Adds H2 boilers as retrofit option for existing gas boilers
   - IMPACT: Missing flexibility option in sector-coupled models
   - DECISION: Could be added as enhancement

5. **Carbon budget constraints** (`set_carbon_constraints()`):
   - Different handling of CO2 budgets for perfect foresight
   - IMPACT: Carbon budget may not work correctly
   - DECISION: Needs investigation

## Validation Checklist

Before using perfect foresight in production:

- [x] Code passes ruff linting
- [x] Python compilation succeeds
- [x] Functions have comprehensive docstrings
- [x] Type hints are complete
- [ ] Integration tests pass (pending test config creation)
- [ ] Results comparison with old workflow (if feasible)
- [ ] Documentation updated in main docs

## Notes for Future Development

1. **Test Configurations**: Create test configs for perfect foresight:
   - `config/test/config.electricity-perfect.yaml`
   - `config/test/config.perfect.yaml`

2. **Phase-Outs**: Consider whether to implement national phase-out policies or handle via config

3. **HVDC Model**: Investigate if HVDC transport model is necessary for perfect foresight line expansion

4. **Heat Pump Efficiency**: Review if current brownfield heat pump efficiency update is sufficient for perfect foresight

5. **Carbon Budgets**: Test carbon budget functionality with perfect foresight

## References

- Original implementation: `scripts/prepare_perfect_foresight.py` (deprecated)
- Reference issue/PR: #1838
- PyPSA documentation: https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#multi-horizon-investment-optimisation
