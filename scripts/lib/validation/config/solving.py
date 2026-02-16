# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Solving configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#solving
"""

from typing import Any

from pydantic import BaseModel, Field, field_validator, model_validator

from scripts.lib.validation.config._base import ConfigModel


class _PostDiscretizationConfig(ConfigModel):
    """Configuration for `solving.options.post_discretization` settings."""

    enable: bool = Field(
        False,
        description="Switch to enable post-discretization of the network. Disabled by default.",
    )
    line_unit_size: float = Field(
        1700, description="Discrete unit size of lines in MW."
    )
    line_threshold: float = Field(
        0.3,
        description="The threshold relative to the discrete line unit size beyond which to round up to the next unit.",
    )
    link_unit_size: dict[str, float] = Field(
        default_factory=lambda: {"DC": 2000, "H2 pipeline": 1200, "gas pipeline": 1500},
        description="Discrete unit size of links in MW by carrier (given in dictionary style).",
    )
    link_threshold: dict[str, float] = Field(
        default_factory=lambda: {"DC": 0.3, "H2 pipeline": 0.3, "gas pipeline": 0.3},
        description="The threshold relative to the discrete link unit size beyond which to round up to the next unit by carrier (given in dictionary style).",
    )
    fractional_last_unit_size: bool = Field(
        False,
        description="When true, links and lines can be built up to p_nom_max. When false, they can only be built up to a multiple of the unit size.",
    )


class _ModelKwargsConfig(BaseModel):
    """Configuration for `solving.options.model_kwargs` settings."""

    solver_dir: str = Field(
        "", description="Absolute path to the directory where linopy saves files."
    )


class _SolvingOptionsConfig(BaseModel):
    """Configuration for `solving.options` settings."""

    clip_p_max_pu: float = Field(
        0.01,
        description="To avoid too small values in the renewables` per-unit availability time series values below this threshold are set to zero.",
    )
    load_shedding: bool | float = Field(
        False,
        description="Add generators with very high marginal cost to simulate load shedding and avoid problem infeasibilities. If load shedding is a float, it denotes the marginal cost in EUR/kWh.",
    )
    curtailment_mode: bool = Field(
        False,
        description="Fixes the dispatch profiles of generators with time-varying p_max_pu by setting `p_min_pu = p_max_pu` and adds an auxiliary curtailment generator (with negative sign to absorb excess power) at every AC bus. This can speed up the solving process as the curtailment decision is aggregated into a single generator per region. Defaults to `false`.",
    )
    noisy_costs: bool = Field(
        True,
        description="Add random noise to marginal cost of generators by :math:`\\mathcal{U}(0.009,0,011)` and capital cost of lines and links by :math:`\\mathcal{U}(0.09,0,11)`.",
    )
    skip_iterations: bool = Field(
        True,
        description="Skip iterating, do not update impedances of branches. Defaults to true.",
    )
    rolling_horizon: bool = Field(
        False,
        description="Switch for rule `solve_operations_network` whether to optimize the network in a rolling horizon manner, where the snapshot range is split into slices of size `horizon` which are solved consecutively. This setting has currently no effect on sector-coupled networks.",
    )
    seed: int = Field(
        123, description="Random seed for increased deterministic behaviour."
    )
    custom_extra_functionality: str | None = Field(
        "../data/custom_extra_functionality.py",
        description="Path to a Python file with custom extra functionality code to be injected into the solving rules of the workflow relative to `rules` directory.",
    )
    io_api: str | None = Field(
        None,
        description="Passed to linopy and determines the API used to communicate with the solver. With the `'lp'` and `'mps'` options linopy passes a file to the solver; with the `'direct'` option (only supported for HIGHS and Gurobi) linopy uses an in-memory python API resulting in better performance.",
    )
    track_iterations: bool = Field(
        False,
        description="Flag whether to store the intermediate branch capacities and objective function values are recorded for each iteration in `network.lines['s_nom_opt_X']` (where `X` labels the iteration)",
    )
    min_iterations: int = Field(
        2,
        description="Minimum number of solving iterations in between which resistance and reactence (`x/r`) are updated for branches according to `s_nom_opt` of the previous run.",
    )
    max_iterations: int = Field(
        3,
        description="Maximum number of solving iterations in between which resistance and reactence (`x/r`) are updated for branches according to `s_nom_opt` of the previous run.",
    )
    transmission_losses: int = Field(
        2,
        description="Add piecewise linear approximation of transmission losses based on n tangents. Defaults to 0, which means losses are ignored.",
    )
    linearized_unit_commitment: bool = Field(
        True,
        description="Whether to optimise using the linearized unit commitment formulation.",
    )
    horizon: int = Field(
        365,
        description="Number of snapshots to consider in each iteration. Defaults to 100.",
    )
    post_discretization: _PostDiscretizationConfig = Field(
        default_factory=_PostDiscretizationConfig,
        description="Post-discretization settings.",
    )
    keep_files: bool = Field(
        False, description="Whether to keep LPs and MPS files after solving."
    )
    store_model: bool = Field(
        False,
        description="Store the linopy model to a NetCDF file after solving. Not supported with rolling_horizon. Not scenario-aware.",
    )
    model_kwargs: _ModelKwargsConfig = Field(
        default_factory=_ModelKwargsConfig, description="Model kwargs for linopy."
    )

    @model_validator(mode="after")
    def check_store_model_rolling_horizon(self):
        if self.rolling_horizon and self.store_model:
            raise ValueError("store_model is not supported with rolling_horizon")
        return self


class _AggPNomLimitsConfig(BaseModel):
    """Configuration for `solving.agg_p_nom_limits` settings."""

    agg_offwind: bool = Field(
        False,
        description="Aggregate together all the types of offwind when writing the constraint (`offwind-all` as a carrier in the `.csv` file). Default is false.",
    )
    agg_solar: bool = Field(
        False,
        description="Aggregate together all the types of electric solar when writing the constraint (`solar-all` as a carrier in the `.csv` file). Default is false.",
    )
    include_existing: bool = Field(
        False,
        description="Take existing capacities into account when writing the constraint. Default is false.",
    )
    file: str = Field(
        "data/agg_p_nom_minmax.csv",
        description="Reference to `.csv` file specifying per carrier generator nominal capacity constraints for individual countries and planning horizons. Defaults to `data/agg_p_nom_minmax.csv`.",
    )


class _ConstraintsConfig(BaseModel):
    """Configuration for `solving.constraints` settings."""

    CCL: bool = Field(
        False,
        description="Add minimum and maximum levels of generator nominal capacity per carrier for individual countries. These can be specified in the file linked at `electricity: agg_p_nom_limits` in the configuration. File defaults to `data/agg_p_nom_minmax.csv`. Does not work with a time resolution resampling.",
    )
    EQ: bool | str = Field(
        False,
        description="Require each country or node to on average produce a minimal share of its total consumption itself. Example: `EQ0.5c` demands each country to produce on average at least 50% of its consumption; `EQ0.5` demands each node to produce on average at least 50% of its consumption.",
    )
    BAU: bool = Field(
        False,
        description="Add a per-`carrier` minimal overall capacity; i.e. at least `40GW` of `OCGT` in Europe; configured in `electricity: BAU_mincapacities`",
    )
    SAFE: bool = Field(
        False,
        description="Add a capacity reserve margin of a certain fraction above the peak demand to which renewable generators and storage do *not* contribute. Ignores network.",
    )


class _SolverConfig(BaseModel):
    """Configuration for `solving.solver` settings."""

    name: str = Field(
        "gurobi",
        description="Solver to use for optimisation problems in the workflow; e.g. clustering and linear optimal power flow.",
    )
    options: str = Field(
        "gurobi-default", description="Link to specific parameter settings."
    )


class _CheckObjectiveConfig(BaseModel):
    """Configuration for `solving.check_objective` settings."""

    enable: bool = Field(False, description="Enable objective value checking.")
    expected_value: float | None = Field(None, description="Expected objective value.")
    atol: float = Field(1_000_000, description="Absolute tolerance.")
    rtol: float = Field(0.01, description="Relative tolerance.")

    @field_validator("expected_value", mode="before")
    @classmethod
    def parse_none_string(cls, v):
        if v == "None" or v is None:
            return None
        return v


class _OETCConfig(BaseModel):
    """Configuration for `solving.oetc` settings (Open Energy Transition Computing cluster support)."""

    name: str = Field(
        "pypsa-eur",
        description="Name identifier for the OETC job.",
    )
    authentication_server_url: str = Field(
        "",
        description="URL of the OETC authentication server for job submission.",
    )
    orchestrator_server_url: str = Field(
        "",
        description="URL of the OETC orchestrator server for job management.",
    )
    cpu_cores: int = Field(
        8,
        description="Number of CPU cores to request for the OETC job. (includes RAM amount at the moment with a factor of 8)",
    )
    disk_space_gb: int = Field(
        50,
        description="Amount of disk space in gigabytes to request for the OETC job.",
    )
    delete_worker_on_error: bool = Field(
        True,
        description="Whether to delete the worker instance when an error occurs during job execution.",
    )


class SolvingConfig(BaseModel):
    """Configuration for `solving` settings."""

    options: _SolvingOptionsConfig = Field(
        default_factory=_SolvingOptionsConfig, description="Solving options."
    )
    agg_p_nom_limits: _AggPNomLimitsConfig = Field(
        default_factory=_AggPNomLimitsConfig,
        description="Aggregate p_nom limits configuration.",
    )
    constraints: _ConstraintsConfig = Field(
        default_factory=_ConstraintsConfig, description="Constraints configuration."
    )
    solver: _SolverConfig = Field(
        default_factory=_SolverConfig, description="Solver configuration."
    )
    solver_options: dict[str, dict[str, Any]] = Field(
        default_factory=lambda: {
            "highs-default": {
                "threads": 1,
                "solver": "ipm",
                "run_crossover": "off",
                "small_matrix_value": 1e-6,
                "large_matrix_value": 1e9,
                "primal_feasibility_tolerance": 1e-5,
                "dual_feasibility_tolerance": 1e-5,
                "ipm_optimality_tolerance": 1e-4,
                "parallel": "on",
                "random_seed": 123,
            },
            "highs-simplex": {
                "solver": "simplex",
                "parallel": "on",
                "primal_feasibility_tolerance": 1e-5,
                "dual_feasibility_tolerance": 1e-5,
                "random_seed": 123,
            },
            "gurobi-default": {
                "threads": 32,
                "method": 2,
                "crossover": 0,
                "BarConvTol": 1e-5,
                "Seed": 123,
                "AggFill": 0,
                "PreDual": 0,
                "GURO_PAR_BARDENSETHRESH": 200,
            },
            "gurobi-numeric-focus": {
                "NumericFocus": 3,
                "method": 2,
                "crossover": 0,
                "BarHomogeneous": 1,
                "BarConvTol": 1e-5,
                "FeasibilityTol": 1e-4,
                "OptimalityTol": 1e-4,
                "ObjScale": -0.5,
                "threads": 8,
                "Seed": 123,
            },
            "gurobi-fallback": {
                "crossover": 0,
                "method": 2,
                "BarHomogeneous": 1,
                "BarConvTol": 1e-5,
                "FeasibilityTol": 1e-5,
                "OptimalityTol": 1e-5,
                "Seed": 123,
                "threads": 8,
            },
            "cplex-default": {
                "threads": 4,
                "lpmethod": 4,
                "solutiontype": 2,
                "barrier.convergetol": 1e-5,
                "feasopt.tolerance": 1e-6,
            },
            "copt-default": {
                "Threads": 8,
                "LpMethod": 2,
                "Crossover": 0,
                "RelGap": 1e-6,
                "Dualize": 0,
            },
            "copt-gpu": {
                "LpMethod": 6,
                "GPUMode": 1,
                "PDLPTol": 1e-5,
                "Crossover": 0,
            },
            "xpress-default": {
                "threads": 8,
                "lpflags": 4,
                "crossover": 0,
                "bargaptarget": 1e-5,
                "baralg": 2,
            },
            "xpress-gpu": {
                "lpflags": 4,
                "crossover": 0,
                "baralg": 4,
                "barhggpu": 1,
                "barhgreltol": 1e-5,
            },
            "cbc-default": {},
            "glpk-default": {},
        },
        description="Dictionaries with solver-specific parameter settings.",
    )

    @field_validator("solver_options")
    @classmethod
    def check_no_gurobi_credentials(cls, v):
        """Prevent Gurobi license credentials from being stored in config."""
        forbidden_keys = {"WLSACCESSID", "WLSSECRET", "LICENSEID"}
        for solver_name, options in v.items():
            if "env" in options:
                found = forbidden_keys & set(options["env"].keys())
                if found:
                    raise ValueError(
                        f"Gurobi license credentials ({', '.join(found)}) must not be set in config to avoid leaking secrets. "
                        "Use a license file instead or check the PyPSA options documentation on how to pass solver_options via environment variables, "
                        'e.g. PYPSA_PARAMS__OPTIMIZE__SOLVER_OPTIONS={"env": {"WLSACCESSID": "...", "WLSSECRET": "...", "LICENSEID": 1234}}'
                    )
        return v

    check_objective: _CheckObjectiveConfig = Field(
        default_factory=_CheckObjectiveConfig,
        description="Objective checking configuration.",
    )
    oetc: _OETCConfig | None = Field(
        None,
        description="Configuration options for Open Energy Transition Computing (OETC) cluster support.",
    )
    mem_mb: int = Field(
        128000,
        description="Estimated maximum memory requirement for solving networks (MB).",
    )
    memory_logging_frequency: int = Field(
        5, description="Interval in seconds at which memory usage is logged."
    )
    runtime: str = Field("48h", description="Runtime in humanfriendly style.")
