# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: CC0-1.0

"""
Test implementation of streamlined PyPSA-EUR workflow compose rules.

This file demonstrates the new 4-step workflow structure:
base → clustered → composed → solved

Based on the implementation plan in STREAMLINE_WORKFLOW_IMPLEMENTATION_PLAN.md
"""


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

    if foresight == "overnight" and len(horizons) > 1:
        raise ValueError(
            "Overnight optimization can only be run for the first horizon."
        )

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


def get_final_horizon():
    """Get the final planning horizon, handling both single values and lists."""
    planning_horizons = config["temporal"]["planning_horizons"]
    if isinstance(planning_horizons, (int, str)):
        return int(planning_horizons)
    else:
        return int(planning_horizons[-1])


def get_solve_inputs(wildcards):
    """Get inputs for solve rule - just the composed network."""
    return {"network": f"networks/{wildcards.run}/composed_{wildcards.horizon}.nc"}


rule base_network:
    output:
        touch("networks/{run}/base.nc"),


rule shapes:
    output:
        touch("resources/{run}/gadm_shapes.geojson"),
        touch("resources/{run}/country_shapes.geojson"),


# Rule to cluster base network (simplified from current cluster_network.py)
rule cluster_network:
    input:
        network="networks/{run}/base.nc",
        gadm_shapes="resources/{run}/gadm_shapes.geojson",
        country_shapes="resources/{run}/country_shapes.geojson",
    output:
        network="networks/{run}/clustered.nc",
        busmap="resources/{run}/busmap.csv",
    params:
        clusters=lambda w: config["clustering"]["cluster_network"]["n_clusters"],
        simplify_network=True,
        focus_weights=None,
    script:
        "../scripts/cluster_network.py"


# Main composition rule - combines all network building steps
rule compose_network:
    input:
        unpack(get_compose_inputs),
    output:
        "networks/{run}/composed_{horizon}.nc",
    params:
        temporal=lambda w: config["temporal"],
        electricity=lambda w: config.get("electricity", {}),
        sector=lambda w: config.get("sector", {}),
        clustering=lambda w: config.get("clustering", {}),
        existing_capacities=lambda w: config.get("existing_capacities", {}),
        baseyear=2023,
        Nyears=5,  # Investment period length
        countries=["DE", "FR", "GB", "IT", "ES", "PL"],
        conventional_carriers=["nuclear", "oil", "OCGT", "CCGT", "coal", "lignite"],
        renewable={"onwind": {"cutout": "europe-2013-era5"}},
        costs={"year": 2030, "version": "v0.6.0", "fill_values": {"FOM": 0}},
        load={"scale": 1.0},
        conventional={"unit_commitment": False},
        time_resolution=None,
        co2limit=None,
        brownfield={"enable": True},
    log:
        "logs/{run}/compose_network_{horizon}.log",
    script:
        "../scripts/compose_network.py"


# Unified solving rule for all foresight approaches
rule solve_network:
    input:
        unpack(get_solve_inputs),
    output:
        "networks/{run}/solved_{horizon}.nc",
    params:
        temporal=lambda w: config["temporal"],
        solver=lambda w: config["solving"]["solver"],
    log:
        "logs/{run}/solve_network_{horizon}.log",
    script:
        "../scripts/solve_network.py"


# Collection rules for different workflow targets
rule prepare_networks:
    """Target rule to prepare all composed networks up to final horizon."""
    input:
        lambda w: expand(
            "networks/{run}/composed_{horizon}.nc",
            run=config["run"]["name"],
            horizon=get_final_horizon(),
        ),


rule solve_networks:
    """Target rule to solve all networks up to final horizon."""
    input:
        lambda w: expand(
            "networks/{run}/solved_{horizon}.nc",
            run=config["run"]["name"],
            horizon=get_final_horizon(),
        ),


rule cluster_networks:
    """Target rule to cluster all base networks."""
    input:
        expand("networks/{run}/clustered.nc", run=config["run"]["name"]),


# Test rules for specific foresight modes
rule test_overnight:
    """Test overnight optimization (single horizon)."""
    input:
        lambda w: f"networks/test-overnight/solved_{get_final_horizon()}.nc",


rule test_myopic:
    """Test myopic optimization (sequential horizons)."""
    input:
        lambda w: f"networks/test-myopic/solved_{get_final_horizon()}.nc",


rule test_perfect:
    """Test perfect foresight optimization (all horizons together)."""
    input:
        lambda w: f"networks/test-perfect/solved_{get_final_horizon()}.nc",


# Validation rule to compare old vs new workflow results
rule validate_results:
    """Compare results between old wildcard-based and new config-driven workflow."""
    input:
        old_network="networks/base_s_128_Co2L-3H_sector.nc",  # Old format
        new_network="networks/{run}/solved_{horizon}.nc",  # New format
    output:
        validation_report="results/{run}/validation_report_{horizon}.html",
    script:
        "../scripts/validate_results.py"


# Cleanup rules
rule clean_networks:
    """Clean all generated network files."""
    shell:
        "rm -rf networks/*/composed_*.nc networks/*/solved_*.nc"


rule clean_test:
    """Clean test-specific files."""
    shell:
        "rm -rf networks/test-*/ logs/test-*/"
