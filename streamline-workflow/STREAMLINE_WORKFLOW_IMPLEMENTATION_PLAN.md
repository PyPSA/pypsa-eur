<!--
SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
SPDX-License-Identifier: CC-BY-4.0
-->
# PyPSA-EUR Workflow Streamlining Implementation Plan

**Reference**: GitHub Discussion #1529  
**Objective**: Refactor PyPSA-EUR from complex wildcard-driven workflow to unified 4-step structure: `base → clustered → composed → solved`

## 📋 Executive Summary

This plan implements the proposal from discussion #1529 to harmonize electricity-only and sector-coupled workflows, eliminate most wildcards in favor of config-driven approach, and reduce workflow complexity while maintaining all functionality.

### Core Design Principles
- **Additive approach**: Keep existing scripts, import functions, combine main sections in new compose script
- **No backward compatibility needed**: Clean break from current approach
- **Config-driven**: Move from wildcards to configuration for file targeting
- **Unified workflows**: Same approach for electricity-only and sector-coupled

## 🎯 Status Update (August 2025)

**✅ PROOF OF CONCEPT COMPLETED**  
**⚠️ PRODUCTION IMPLEMENTATION IN PROGRESS WITH ISSUES**

The streamlined workflow concept has been successfully tested via `rules/test_streamlined.smk`, but the production implementation in `rules/compose.smk` currently has syntax errors that need resolution.

### Current Implementation Status
| Component | Status | Details |
|-----------|--------|---------|
| **4-Step Workflow** | ✅ Tested | `base → clustered → composed_{horizon} → solved_{horizon}` |
| **Config-Driven Approach** | ✅ Tested | No wildcard expansion, `temporal` section replaces `scenario` |
| **Multi-Foresight Support** | ✅ Tested | Overnight, Myopic, Perfect foresight modes validated |
| **Test Implementation** | ✅ Working | `rules/test_streamlined.smk` functional |
| **Production Rules** | ⚠️ Broken | `rules/compose.smk` has syntax errors (line 54-57, 161) |
| **Compose Script** | ✅ Exists | `scripts/compose_network.py` (23KB) implemented |
| **Test Configs** | ✅ Complete | Both original and streamlined configs exist |

### Recent Commits (streamline-workflow branch)
```
279b017e refac: brute force rename paths, remove wildcards, consolidate solve rule
7f2ef1a9 add streamlined configs  
051060db Implement streamlined workflow with compose_network.py
685e1684 Refine implementation plan based on successful testing
16bdb75d Implement streamlined workflow test framework
```

### Existing Files and Artifacts
#### Test Implementation (Working)
- **`rules/test_streamlined.smk`**: Complete test workflow (204 lines)
- **`config/test/config.overnight.yaml`**: Original test config
- **`config/test/config.myopic.yaml`**: Original test config
- **`config/test/config.perfect.yaml`**: Original test config

#### Production Implementation (In Progress)  
- **`rules/compose.smk`**: Production rules (10KB) - **HAS SYNTAX ERRORS**
- **`scripts/compose_network.py`**: Main composition script (23KB) - Implemented
- **`config/test/config.streamlined-overnight.yaml`**: Streamlined config
- **`config/test/config.streamlined-myopic.yaml`**: Streamlined config  
- **`config/test/config.streamlined-perfect.yaml`**: Streamlined config

### Critical Issues Fixed (August 16, 2025)
1. **✅ FIXED: Syntax Error in `compose.smk`**: Dict construction mixing positional/keyword args
   - Restructured to use `inputs.update()` pattern for proper merging
2. **✅ FIXED: Duplicate Rule Names**: 
   - `build_population_layouts` → `build_clustered_population_layouts`
   - `build_co2_sequestration_potentials` → `build_clustered_co2_sequestration_potentials`
3. **✅ FIXED: Missing Rule File**: Changed `solve_electricity.smk` → `solve.smk` in Snakefile
4. **✅ FIXED: Wildcard Consistency**: Fixed log/benchmark wildcards in postprocess.smk rules

### Key Learnings from Testing
1. **Horizon handling needs flexibility**: Single values vs lists require helper functions
2. **Config robustness essential**: Missing sections cause workflow failures 
3. **Sequential dependencies work**: Myopic and perfect foresight DAGs validated
4. **Rule structure scales**: Same rules work across all foresight modes
5. **Lambda functions required**: Dynamic horizon resolution needs runtime evaluation

---

## 🎯 Target Architecture

### New 4-Step Workflow
```
1. base       → networks/{run}/base.nc
2. clustered  → networks/{run}/clustered.nc  
3. composed   → networks/{run}/composed_{horizon}.nc
4. solved     → networks/{run}/solved_{horizon}.nc (last horizon = final result)
```

### Current vs Target File Naming
```
CURRENT:  networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc
TARGET:   networks/{run}/solved_{last_horizon}.nc
```

All configuration (clusters, opts, sector_opts, horizons) moves from wildcards to scenario config.

---

## 📊 Current State Analysis

### Wildcard Usage Mapping
- `{clusters}`: 37, 128, 256, etc. → `config.clustering.cluster_network.n_clusters`
- `{opts}`: electricity constraints → `config.electricity.extendable_carriers`, etc.
- `{sector_opts}`: sector settings → `config.sector.*`  
- `{planning_horizons}`: 2030, 2050 → `config.temporal.planning_horizons`

### Key Scripts Dependencies
| Script | Main Functions | Snakemake Dependencies |
|--------|---------------|----------------------|
| `add_electricity.py` | `attach_load()`, `attach_wind_and_solar()`, `attach_storageunits()` | base_network, tech_costs, powerplants |
| `prepare_sector_network.py` | `add_heat()`, `add_transport()`, `add_industry()` | elec network, industrial demand, heat demand |
| `add_existing_baseyear.py` | `add_power_capacities_installed_before_baseyear()` | network, powerplants, heating data |
| `prepare_network.py` | `average_every_nhours()`, `add_co2limit()` | clustered network, tech costs |
| `cluster_network.py` | `busmap_for_n_clusters()`, `clustering_for_n_clusters()` | simplified network, admin shapes |
| `simplify_network.py` | `simplify_network_to_380()`, `remove_stubs()` | base network, admin shapes |

---

## 🔧 Implementation Strategy

### Phase 1: Infrastructure Setup

#### 1.1 Create New Scripts
```bash
scripts/compose_network.py     # Main composition orchestrator
rules/compose.smk             # New workflow rules
# Modify existing scripts/solve_network.py for unified approach
```

#### 1.2 Configuration Schema Updates
```yaml
# Remove scenario section entirely (no more wildcards!)
# scenario:  # REMOVED - no more wildcard expansion

# Add new temporal section for planning settings
temporal:
  planning_horizons: [2030, 2040, 2050]
  foresight: "myopic"  # or "perfect", "overnight"

# Everything else uses existing config sections directly:
# - clustering.cluster_network.n_clusters (was {clusters})
# - electricity.* settings (was {opts})
# - sector.* settings (was {sector_opts})
```

#### 1.2.1 Validated Configuration Pattern (from Testing)
Based on our testing, the following configuration structure has been validated:

```yaml
# Required: Run identification for file paths
run:
  name: "test-sector-overnight"  # Maps to {run} wildcard

# Required: Temporal planning configuration  
temporal:
  foresight: overnight  # or "myopic", "perfect"
  planning_horizons: 2040  # Single value OR [2030, 2040, 2050] list

# Required: Clustering configuration with explicit n_clusters
clustering:
  temporal:
    resolution_sector: 24h
  cluster_network:
    algorithm: kmeans
    n_clusters: 5  # Small number for testing
    hac_features:
    - wnd100m
    - influx_direct

# Existing sections work as-is (electricity, sector, solving, etc.)
```

**Key Testing Insights:**
- `planning_horizons` must handle both single values and lists flexibly
- `cluster_network.n_clusters` is required (was previously derived from wildcard)
- Missing config sections need graceful defaults using `.get()` methods
- `run.name` provides the directory structure for organized outputs

#### 1.3 New Workflow Rules Structure
```python
# rules/compose.smk
rule cluster_network:
    input: "networks/{run}/base.nc"
    output: "networks/{run}/clustered.nc"
    script: "../scripts/cluster_network.py"  # Use existing script with modifications

rule compose_network:
    input: unpack(get_compose_inputs)  # Dynamic inputs based on horizon
    output: "networks/{run}/composed_{horizon}.nc"
    script: "../scripts/compose_network.py"

rule solve_network:
    input: "networks/{run}/composed_{horizon}.nc"
    output: "networks/{run}/solved_{horizon}.nc"
    script: "../scripts/solve_network.py"
```

#### 1.4 Validated Implementation Patterns (from Testing)

**Dynamic Input Resolution:**
```python
def get_compose_inputs(wildcards):
    """Determine inputs for compose rule based on foresight and horizon."""
    temporal = config["temporal"]
    foresight = temporal["foresight"]
    horizon = int(wildcards.horizon)
    
    # Handle both single value and list for planning_horizons
    planning_horizons = temporal["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        horizons = [int(planning_horizons)]
    else:
        horizons = [int(h) for h in planning_horizons]

    inputs = {"clustered": f"networks/{wildcards.run}/clustered.nc"}

    if horizon != horizons[0]:
        # Not first horizon - need previous network
        prev_horizon = horizons[horizons.index(horizon) - 1]

        if foresight == "myopic":
            # Myopic uses solved network from previous horizon
            inputs["network_previous"] = (
                f"networks/{wildcards.run}/solved_{prev_horizon}.nc"
            )
        else:  # perfect foresight
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = (
                f"networks/{wildcards.run}/composed_{prev_horizon}.nc"
            )

    return inputs
```

**Strict Parameter Handling:**
```python
# In rule params use the config_handler function 
params:
    temporal=config_handler("temporal"),
    electricity=config_handler("electricity"),
    sector=config_handler("sector"),
    clustering=config_handler("clustering"),
    another_key=config_handler("electricity", "another_key"),
```

**Dynamic Horizon Resolution:**
```python
def get_final_horizon():
    """Get the final planning horizon, handling both single values and lists."""
    planning_horizons = config["temporal"]["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        return int(planning_horizons)
    else:
        return int(planning_horizons[-1])

# Use in lambda functions for rule inputs
rule test_overnight:
    input:
        lambda w: f"networks/test-overnight/solved_{get_final_horizon()}.nc",
```

### Phase 2: Compose Network Implementation

#### 2.1 `scripts/compose_network.py` Structure
```python
"""
Compose network by importing functions from existing scripts
and combining their main section logic directly.
"""

# Import all required functions from existing scripts
from add_electricity import (
    attach_load, attach_wind_and_solar, attach_storageunits,
    attach_conventional_generators, sanitize_carriers,
    load_costs, load_and_aggregate_powerplants
)
from prepare_sector_network import (
    add_heat, add_transport, add_industry, add_co2_tracking,
    patch_electricity_network, define_spatial, add_generation,
    add_storage_and_grids, cluster_heat_buses
)
from add_existing_baseyear import (
    add_build_year_to_new_assets, 
    add_power_capacities_installed_before_baseyear,
    add_heating_capacities_installed_before_baseyear
)
from add_brownfield import add_brownfield
from prepare_network import (
    average_every_nhours, add_co2limit, add_emission_prices,
    set_transmission_limit, apply_time_segmentation
)

def main():
    """Main composition logic - directly using imported functions as in original scripts."""
    # Load configuration and inputs
    config = snakemake.config
    params = snakemake.params
    
    # Determine current planning horizon and foresight mode
    current_horizon = int(snakemake.wildcards.horizon)
    temporal_config = config["temporal"]
    foresight = temporal_config["foresight"]
    horizons = temporal_config["planning_horizons"]
    
    # Load base network - either clustered.nc or previous horizon's output
    if current_horizon != horizons[0]:
        # Multi-horizon case: load previous network
        if foresight == "myopic":
            # For myopic, input is solved network from previous horizon
            network = pypsa.Network(snakemake.input.network_previous)
            network_current = pypsa.Network(snakemake.input.clustered)
            
            # Apply brownfield constraints from previous solved network
            add_brownfield(network, network_current, params.brownfield)
            network = network_current
        else:  # perfect foresight
            # For perfect, input is composed network from previous horizon
            network_previous = pypsa.Network(snakemake.input.network_previous)
            network_current = pypsa.Network(snakemake.input.clustered)
            
            # Concatenate networks (combine previous horizons with current)
            network = concatenate_networks(network_previous, network_current, current_horizon)
    else:
        # First horizon or overnight - start from clustered network
        network = pypsa.Network(snakemake.input.clustered)
    
    # Step 1: Add electricity components (from add_electricity.py main section)
    Nyears = params.Nyears
    costs = load_costs(snakemake.input.tech_costs, params.costs, Nyears)
    
    ppl, ppl_map = load_and_aggregate_powerplants(
        snakemake.input.powerplants,
        params.countries,
        params.conventional_carriers
    )
    
    attach_load(
        network,
        snakemake.input.load,
        snakemake.input.busmap,
        params.countries,
        params.load
    )
    
    attach_conventional_generators(
        network,
        snakemake.input,
        costs,
        params.conventional
    )
    
    attach_wind_and_solar(
        network,
        snakemake.input,
        costs,
        params.renewable
    )
    
    if params.electricity["extendable_carriers"]:
        attach_storageunits(network, costs, params.electricity)
    
    # Step 2: Add sector components if enabled (from prepare_sector_network.py main section)
    if params.sector["enabled"]:
        patch_electricity_network(network, params)
        
        spatial = define_spatial(network.buses.index, params)
        
        add_co2_tracking(network, params.sector)
        
        add_generation(network, costs, params)
        
        add_storage_and_grids(network, costs, spatial, params)
        
        add_heat(network, snakemake.input, costs, params)
        
        add_transport(network, snakemake.input, costs, params) 
        
        add_industry(network, snakemake.input, costs, params)
        
        if params.sector["cluster_heat_buses"]:
            cluster_heat_buses(network, params.sector)
    
    # Step 3: Add existing capacities if baseyear (from add_existing_baseyear.py)
    if params.existing_capacities["enabled"]:
        add_build_year_to_new_assets(network, params.baseyear)
        
        add_power_capacities_installed_before_baseyear(
            network, 
            params.baseyear,
            snakemake.input.powerplants,
            params.existing_capacities
        )
        
        add_heating_capacities_installed_before_baseyear(
            network,
            params.baseyear,
            snakemake.input,
            params.existing_capacities
        )
    
    # Step 4: Apply temporal resolution and constraints (from prepare_network.py)
    if params.time_resolution:
        network = average_every_nhours(network, params.time_resolution)
    
    if params.co2limit:
        add_co2limit(network, params.co2limit, Nyears)
        
    if params.electricity["transmission_limit"]:
        set_transmission_limit(network, params.electricity["transmission_limit"])
    
    add_emission_prices(network, snakemake.input.co2_price, params)
    
    # Step 5: Final processing
    sanitize_carriers(network, config)
    
    # Export composed network
    network.export_to_netcdf(snakemake.output[0])

def concatenate_networks(network_previous, network_current, current_horizon):
    """Concatenate previous horizons' network with current horizon."""
    # For perfect foresight: combine all components from all horizons
    # Previous network contains horizons [2030, 2035, ...]
    # Current network contains only current horizon [2040]
    # Result should contain [2030, 2035, 2040]
    
    # Create new network with multi-indexed snapshots
    import pandas as pd
    
    # Get all planning horizons from previous network
    previous_horizons = network_previous.meta.get("horizon", [])
    
    # Add current horizon
    all_horizons = previous_horizons + [current_horizon]
    
    # Create multi-period network structure
    # This involves:
    # 1. Combining snapshots with horizon index
    # 2. Merging components (buses, generators, etc.) with horizon suffixes
    # 3. Handling inter-temporal constraints and variables
    
    # Implementation depends on PyPSA multi-period API
    return network  # TODO: Complete implementation
```

#### 2.2 Multi-Horizon DAG Design

Based on FabianHofmann's comment, implement different workflows:

**Overnight Workflow**:
```
compose_{horizon} → solve_{horizon} → solved.nc
```
- Single horizon (e.g., 2050)
- Direct optimization of composed network
- Final output is copy of single horizon result

**Myopic Workflow**:
```
compose_2030 (input: clustered) → solve_2030 → compose_2040 (input: solved_2030) → solve_2040 → ...
```
- Each horizon is composed with brownfield constraints from previous solved network
- Each horizon is solved independently
- Final output contains last horizon's solved network

**Perfect Foresight Workflow**:
```
compose_2030 → compose_2040 → ... → compose_2050 → solve_2050 (all horizons)
```
- Each horizon is composed and concatenated with previous composed networks
- Only the last horizon's solve rule runs (solves all horizons together)
- Earlier horizons are never individually solved (not needed since solve_X files aren't requested)
- Final output contains all horizons in multi-period network

**Key Differences**:
| Aspect | Myopic | Perfect Foresight |
|--------|--------|------------------|
| Input to compose | Previous solved network | Previous composed network |
| Brownfield logic | Applied in compose step | Not needed (all horizons optimized together) |
| Network structure | Single horizon per file | Multi-horizon concatenated |
| Solving | Each horizon separately (simple optimization) | All horizons together |
| Output | Sequential solved networks | Single multi-period network |

#### 2.3 `rules/compose.smk` Unified Rules with Horizon Wildcard
```python
def get_compose_inputs(wildcards):
    """Determine inputs for compose rule based on foresight and horizon."""
    config = get_config()
    temporal = config["temporal"]
    foresight = temporal["foresight"]
    horizon = int(wildcards.horizon)
    horizons = temporal["planning_horizons"]
    
    inputs = {"clustered": f"networks/{wildcards.run}/clustered.nc"}
    
    if horizon != horizons[0]:
        # Not first horizon - need previous network
        prev_horizon = horizons[horizons.index(horizon) - 1]
        
        if foresight == "myopic":
            # Myopic uses solved network from previous horizon
            inputs["network_previous"] = f"networks/{wildcards.run}/solved_{prev_horizon}.nc"
        else:  # perfect foresight
            # Perfect foresight uses composed network from previous horizon
            inputs["network_previous"] = f"networks/{wildcards.run}/composed_{prev_horizon}.nc"
    
    # Add other required inputs using resource functions with {run} paths
    inputs.update({
        "tech_costs": resources("costs_{horizon}.csv"),
        "powerplants": resources("powerplants.csv"),
        "load": resources("electricity_demand.csv"),
        "busmap": resources("busmap.csv"),
        "profile_solar": resources("profile_{run}/solar.nc"),
        "profile_onwind": resources("profile_{run}/onwind.nc"),
        # ... other inputs as needed, some with {horizon} wildcard
    })
    
    return inputs

rule compose_network:
    input: unpack(get_compose_inputs)
    output: "networks/{run}/composed_{horizon}.nc"
    params:
        temporal=config["temporal"],
        electricity=config["electricity"],
        sector=config["sector"],
        clustering=config["clustering"],
        existing_capacities=config["existing_capacities"],
        # Cherry-pick other config sections as needed
    script: "../scripts/compose_network.py"

rule solve_network:
    input: "networks/{run}/composed_{horizon}.nc"
    output: "networks/{run}/solved_{horizon}.nc"
    params:
        temporal=config["temporal"],
        solver=config["solving"]["solver"]
    script: "../scripts/solve_network.py"

# No solve_final rule needed! 
# For perfect foresight: solve_{last_horizon} takes composed_{last_horizon} (all horizons) and solves once
# For myopic/overnight: solve_{last_horizon} is the natural final result
```

### Phase 3: Solving Integration

#### 3.1 Simplified Solving Script
```python
# scripts/solve_network.py (modify existing script)
"""Simplified unified solving for all foresight approaches."""

def main():
    """Main solving logic - simplified for all cases."""
    # Load network
    network = pypsa.Network(snakemake.input[0])
    
    # Get config and params
    config = snakemake.config
    params = snakemake.params
    
    # Prepare solver options
    solver_options = params.solver
    solver_name = solver_options["solver_name"]
    
    # Determine if we need to pass snapshots (for myopic only)
    temporal = params.temporal
    foresight = temporal["foresight"]
    current_horizon = int(snakemake.wildcards.horizon)
    
    if foresight == "myopic":
        # For myopic: solve only current horizon's snapshots
        snapshots = network.snapshots[network.snapshots.get_level_values("period") == current_horizon]
        status, condition = network.optimize(
            solver_name=solver_name,
            solver_options=solver_options,
            snapshots=snapshots
        )
    else:
        # For overnight and perfect foresight: solve all snapshots
        status, condition = network.optimize(
            solver_name=solver_name,
            solver_options=solver_options
        )
    
    # Check optimization status
    if status != "ok":
        logger.warning(f"Optimization failed with status {status}")
    
    # Export solved network
    network.export_to_netcdf(snakemake.output[0])

if __name__ == "__main__":
    main()
```

### Phase 4: Integration and Testing

#### 4.1 Rule Dependencies Update
```python
# Update rules/collect.smk with new targets
rule solve_networks:
    input:
        lambda w: expand(
            "networks/{run}/solved_{horizon}.nc",
            run=config["run"]["name"],  # Note: run.name can be a list
            horizon=config["temporal"]["planning_horizons"][-1]
        )

rule prepare_networks:
    input: 
        lambda w: expand(
            "networks/{run}/composed_{horizon}.nc",
            run=config["run"]["name"],  # Note: run.name can be a list
            horizon=config["temporal"]["planning_horizons"][-1]
        )

rule cluster_networks:
    input:
        expand(
            "networks/{run}/clustered.nc",
            run=config["run"]["name"]  # Note: run.name can be a list
        )
```

#### 4.2 Configuration Migration
- Remove `scenario` section entirely (no wildcard expansion)
- Add new `temporal` section for planning horizons and foresight
- Use existing config sections directly (clustering, electricity, sector)
- No migration script needed - clean break from wildcards

#### 4.3 Testing Strategy
```bash
# Basic workflow validation
snakemake cluster_networks -n  # Dry run clustering step
snakemake prepare_networks -n  # Dry run compose step
snakemake solve_networks -n    # Dry run solve step

# Test each foresight mode specifically
snakemake networks/test-run/solved_2050.nc -n  # Overnight mode
snakemake networks/test-myopic/solved_2050.nc -n  # Myopic mode  
snakemake networks/test-perfect/solved_2050.nc -n  # Perfect foresight mode

# Validate results match current workflow  
python scripts/validate_results.py \
  --old networks/base_s_128_Co2L-3H_sector.nc \
  --new networks/{run}/solved_{last_horizon}.nc

# Test horizon-specific input handling
snakemake networks/test-run/composed_2030.nc -n  # First horizon
snakemake networks/test-run/composed_2040.nc -n  # Middle horizon
snakemake networks/test-run/composed_2050.nc -n  # Last horizon
```

---

## 🔄 Implementation Timeline

### ✅ Phase 0: Proof of Concept (COMPLETED - July 2025)
- [x] **Create test workflow framework** (`rules/test_streamlined.smk`)
- [x] **Validate 4-step architecture** with all foresight modes
- [x] **Test configuration schema** with `temporal` section replacing `scenario`
- [x] **Implement dynamic dependencies** for myopic/perfect foresight 
- [x] **Robust parameter handling** with defaults for missing config sections
- [x] **Multi-horizon support** with flexible planning_horizons handling
- [x] **End-to-end workflow validation** via dry-run testing

### ⚠️ Phase 1: Script Implementation (CURRENT - NEEDS FIXES)
- [x] ~~Create `scripts/compose_network.py` based on validated test framework~~ **DONE**
- [x] ~~Import and integrate functions from existing scripts~~ **DONE** 
- [ ] **FIX CRITICAL**: Resolve syntax errors in `rules/compose.smk`
- [ ] Implement network concatenation for perfect foresight
- [ ] Add brownfield logic for myopic optimization  
- [ ] Test with real input data and actual script execution
- [ ] Verify integration with existing rule dependencies

### Phase 2: Integration & Production (BLOCKED BY PHASE 1)
- [ ] Update main `Snakefile` to include compose rules
- [ ] Create collection rules (`solve_networks`, `prepare_networks`)
- [ ] Migrate configuration templates and examples
- [ ] Comprehensive testing with full scenarios
- [ ] Performance optimization and memory management
- [ ] Resolve conflicts with existing workflow rules

### Phase 3: Documentation & Migration
- [ ] Update user documentation with new workflow approach
- [ ] Create migration guide for existing users
- [ ] Add configuration validation and helpful error messages
- [ ] Community feedback and iteration

## 🚨 Immediate Actions Required

### 1. Fix Syntax Errors in `compose.smk`
The dict construction in `get_compose_inputs()` has mixed positional and keyword arguments:
```python
# BROKEN: Lines 54-57
inputs = dict(
    unpack(input_profile_tech),  # OK - positional at start
    base_network=resources("..."),  # Keyword args
    # ...
    unpack(input_profile_offwind),  # ERROR: positional after keyword
    unpack(input_heat_source_power),  # ERROR: positional after keyword
    **rules.cluster_gas_network.output,  # Should use proper syntax
)
```

**Solution**: Move all unpacked/positional arguments to the beginning or use `.update()` method.

### 2. Verify Rule Dependencies
- Check if `cluster_gas_network` and `build_gas_input_locations` rules exist
- Ensure all referenced input functions (`input_profile_tech`, etc.) are defined
- Validate that `config_provider` and `resources` functions work correctly

### 3. Test Workflow Execution
After fixing syntax errors:
```bash
# Test compose workflow with streamlined configs
snakemake -s rules/compose.smk networks/test-run/composed_2040.nc \
  --configfile config/test/config.streamlined-overnight.yaml -n

# Verify all foresight modes
for mode in overnight myopic perfect; do
  snakemake -s rules/compose.smk test_${mode} \
    --configfile config/test/config.streamlined-${mode}.yaml -n
done
```

### 4. Integration Considerations
- The workflow appears to have moved beyond the original plan with additional complexity
- Recent commits show "brute force rename paths" and "remove wildcards" indicating aggressive refactoring
- Need to ensure backward compatibility or clear migration path

---

## 🎯 Success Criteria

### ✅ Proof of Concept Achievements (COMPLETED)
1. **✅ Simplified Structure**: 4-step workflow with only `{run}` and `{horizon}` wildcards validated
2. **✅ Config-Driven**: Target files determined by configuration, not wildcards - working
3. **✅ Unified DAG**: Single rule structure works for overnight, myopic, and perfect foresight - tested
4. **✅ Robust Parameter Handling**: Graceful defaults for missing config sections - implemented
5. **✅ Multi-Foresight Support**: All three foresight modes working with correct dependencies

### 🎯 Production Implementation Goals (NEXT)
1. **Functional Equivalence**: New workflow produces identical results to current workflow
2. **Unified Approach**: Same workflow handles electricity-only and sector-coupled
3. **Maintainable**: Clear separation of concerns, reusable functions
4. **Performance**: Efficient memory usage and execution time
5. **User-Friendly**: Clear error messages and configuration validation

---

## 🚨 Risk Mitigation

### Technical Risks
- **Large intermediate files**: Monitor disk usage, implement cleanup rules
- **Memory usage**: Test with largest scenarios, optimize data loading
- **Function compatibility**: Extensive unit testing of imported functions

### Process Risks  
- **Regression**: Comprehensive result validation pipeline
- **Timeline**: Start with minimal viable implementation, iterate
- **Coordination**: Regular check-ins, clear interfaces between components

---

## 🔍 Open Ends and Potential Issues

### 🚩 Current State Assessment
The implementation has diverged from the original plan with more aggressive changes than initially proposed:

1. **Scope Creep**: Recent commits show "brute force rename paths" suggesting more invasive changes than the additive approach originally planned
2. **Broken State**: Production rules have syntax errors preventing execution
3. **Incomplete Testing**: While test framework exists, actual execution with real data hasn't been validated
4. **Missing Integration**: No clear path from current broken state to working production

### ⚠️ Technical Debt and Risks

#### Immediate Issues
- **Syntax Errors**: `compose.smk` cannot execute due to dict construction errors
- **Undefined References**: Rules reference `cluster_gas_network` and `build_gas_input_locations` - need verification
- **Path Inconsistencies**: Mix of old wildcard paths and new streamlined paths

#### Integration Challenges  
- **Parallel Development**: Main branch continues evolving while streamline-workflow diverges
- **Testing Gap**: No evidence of end-to-end execution with actual data
- **Config Migration**: Two sets of test configs (original vs streamlined) without clear relationship

#### Architecture Concerns
- **Aggressive Refactoring**: "Brute force" approach may have broken assumptions
- **Missing Backwards Compatibility**: No migration path for existing users
- **Unclear Dependencies**: Rule dependencies between compose and existing rules unclear

## 🔍 Open Questions for Discussion

### ✅ Resolved Questions (from Testing)
1. **~~Rule Dependencies~~**: ✅ **SOLVED** - Dynamic input functions with lambda expressions work perfectly
2. **~~Configuration Flexibility~~**: ✅ **SOLVED** - Helper functions handle single values vs lists robustly  
3. **~~Foresight Mode Implementation~~**: ✅ **SOLVED** - Single rule structure works for all modes
4. **~~Parameter Robustness~~**: ✅ **SOLVED** - `.get()` methods with defaults prevent config errors

### 🤔 Critical Unresolved Questions

1. **Network Concatenation Implementation**: What's the best approach for concatenating networks across horizons in PyPSA?
   - How to handle multi-indexed snapshots for perfect foresight?
   - How to merge components with horizon-specific attributes?
   - Should we use PyPSA's native multi-period capabilities or custom implementation?

2. **Input Data Handling**: How should time-varying input data be handled across horizons?
   - Renewable profiles for different years (weather data variability)
   - Demand projections and growth rates
   - Cost trajectories and technology learning curves
   - Should inputs be horizon-specific in file paths?

3. **Brownfield Implementation Details**: For myopic optimization:
   - Should all existing capacities be locked or allow decommissioning?
   - How to handle technology-specific lifetime constraints?
   - Interface between `add_brownfield` and composed network structure
   - How to transfer capacity constraints between horizons?

4. **Memory Management**: For perfect foresight with many horizons:
   - Expected memory usage for concatenated networks (3+ horizons)
   - Strategies for handling very large multi-period networks
   - Should we implement checkpointing between horizons?
   - Disk vs memory tradeoffs for intermediate files

5. **Production Integration**: How to integrate with existing workflow?
   - Should we add to main `Snakefile` or keep separate?
   - How to handle backward compatibility during transition?
   - Configuration migration strategy for existing users

6. **Backward Compatibility Bridge**: Although not required, should we provide:
   - Script to convert old wildcard-based configs to new format?
   - Mapping table for users to understand the transition?

---

## 📚 Recommended Recovery Path

### Option 1: Fix and Continue (Recommended)

1. **Fix Syntax Errors**: Resolve dict construction issues in `compose.smk`
2. **Test Incrementally**: Start with the simplest case (overnight) before complex (myopic/perfect)
3. **Verify Dependencies**: Ensure all referenced rules and functions exist
4. **Document Changes**: Update this plan with actual implementation details
5. **Create Tests**: Add unit tests for compose_network.py functions

### Option 2: Revert to Test Implementation

1. **Use Working Test**: Start from `test_streamlined.smk` which is known to work
2. **Gradual Migration**: Move functionality piece by piece to production
3. **Maintain Stability**: Keep test and production separate until fully validated
4. **Incremental Integration**: Merge with main workflow only after full validation

### Option 3: Re-evaluate Approach

1. **Review Original Goals**: Ensure current implementation aligns with discussion #1529
2. **Consider Alternatives**: Less invasive approaches that maintain compatibility
3. **Phased Rollout**: Implement as optional alternative workflow first
4. **Community Input**: Get feedback before committing to breaking changes

## 🎯 Priority Actions (Updated August 16, 2025)

### ✅ Completed Actions
1. **Fixed syntax errors** in compose.smk - dict construction restructured
2. **Fixed duplicate rule names** in build_sector.smk - renamed to clustered variants
3. **Fixed missing rule file** - corrected solve_electricity.smk to solve.smk
4. **Fixed wildcard consistency** in postprocess.smk

### 🚧 Remaining Actions
1. **SHORT-TERM**: Complete network concatenation for perfect foresight
2. **SHORT-TERM**: Implement brownfield logic for myopic optimization  
3. **MEDIUM-TERM**: Validate with real data (not just test configs)
4. **LONG-TERM**: Full migration guide and documentation

## 📊 Success Metrics

- [ ] Compose workflow executes without errors
- [ ] Test cases pass for all three foresight modes  
- [ ] Results match current workflow output
- [ ] Performance is comparable or better
- [ ] Clear migration path documented

This plan provides a comprehensive assessment of the current state and a roadmap for completing the workflow streamlining implementation.

