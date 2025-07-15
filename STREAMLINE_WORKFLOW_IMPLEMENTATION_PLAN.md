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

### Week 1-2: Infrastructure
- [ ] Create `scripts/compose_network.py` with basic structure
- [ ] Create `rules/compose.smk` with single-horizon rules  
- [ ] Set up configuration schema extensions
- [ ] Test basic composition for electricity-only case

### Week 3-4: Multi-Horizon Logic
- [ ] Implement myopic workflow DAG
- [ ] Implement perfect foresight workflow DAG
- [ ] Add horizon merging and brownfield logic
- [ ] Test multi-horizon composition

### Week 5-6: Solving Integration
- [ ] Create `scripts/solve_unified.py`
- [ ] Integrate all three foresight approaches
- [ ] Update solve rules in workflow
- [ ] End-to-end testing

### Week 7-8: Integration & Cleanup
- [ ] Update collection rules and targets
- [ ] Create configuration migration tools
- [ ] Comprehensive testing across scenarios
- [ ] Documentation updates

---

## 🎯 Success Criteria

1. **Functional Equivalence**: New workflow produces identical results to current workflow
2. **Simplified Structure**: 4-step workflow with only `{run}` and `{horizon}` wildcards
3. **Unified Approach**: Same workflow handles electricity-only and sector-coupled
4. **Config-Driven**: Target files determined by configuration, not wildcards
5. **Maintainable**: Clear separation of concerns, reusable functions
6. **Unified DAG**: Single rule structure works for overnight, myopic, and perfect foresight

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

## 🔍 Open Questions for Discussion

1. **Network Concatenation Implementation**: What's the best approach for concatenating networks across horizons in PyPSA?
   - How to handle multi-indexed snapshots?
   - How to merge components with horizon-specific attributes?
   - Should we use PyPSA's native multi-period capabilities or custom implementation?

2. **Input Data Handling**: How should time-varying input data be handled across horizons?
   - Renewable profiles for different years
   - Demand projections
   - Cost trajectories

3. **Brownfield Implementation Details**: For myopic optimization:
   - Should all existing capacities be locked or allow decommissioning?
   - How to handle technology-specific lifetime constraints?
   - Interface between `add_brownfield` and composed network structure

4. **Memory Management**: For perfect foresight with many horizons:
   - Expected memory usage for concatenated networks
   - Strategies for handling very large multi-period networks
   - Should we implement checkpointing between horizons?

5. **Rule Dependencies**: How to elegantly express the DAG dependencies?
   - Dynamic rule generation based on foresight mode?
   - Separate .smk files for myopic vs perfect foresight?
   - How to handle the final aggregation step?

6. **Backward Compatibility Bridge**: Although not required, should we provide:
   - Script to convert old wildcard-based configs to new format?
   - Mapping table for users to understand the transition?

---

## 📚 Next Steps

1. **Review and Approve Plan**: Get team feedback on approach and timeline
2. **Start Implementation**: Begin with Phase 1 infrastructure setup  
3. **Iterative Development**: Regular testing and refinement
4. **Documentation**: Update user guides and developer documentation
5. **Community Communication**: Share progress and gather feedback

This plan provides a comprehensive roadmap for implementing the workflow streamlining while maintaining all current functionality and providing a foundation for future modularity improvements.