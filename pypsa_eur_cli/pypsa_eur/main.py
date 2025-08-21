import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import typer
import yaml
from pypsa_eur.cli_helpers import (
    check_in_project_root,
    comma_separated_list,
    generate_job_id,
    recursive_merge,
)
from pypsa_eur.progress_tracker import ProgressTracker
from snakemake.utils import validate
from tqdm import tqdm

app = typer.Typer(
    help="PyPSA-Eur CLI: Command Line Interface for European Power System Analysis",
    no_args_is_help=True,
)


@app.command(name="run")
def run(
        name: str = typer.Option(
            None,
            "--name",
            "-n",
            help="Name for this run (will be used as job ID).",
        ),
        configfile: Path = typer.Option(
            None,
            "--configfile",
            "-c",
            help="Path to configuration file (default: config/config.default.yaml).",
        ),
        verbose: bool = typer.Option(
            False,
            "--verbose",
            "-v",
            help="Show live Snakemake output during workflow execution.",
        ),
        cores: int = typer.Option(
            None,
            "--cores",
            help="Number of cores to use for Snakemake. If not specified, all available cores will be used.",
        ),
        snapshots: str = typer.Option(
            None,
            "--snapshots",
            "-s",
            help="Override snapshots configuration (format: 'start,end' e.g., '2013-01-01,2014-01-01').",
        ),
        planning_horizons: str = typer.Option(
            None,
            "--planning-horizons",
            help="Override planning horizons as comma-separated list (e.g., 2030,2040,2050).",
        ),
        clusters: str = typer.Option(
            None,
            "--clusters",
            help="Override network clusters as comma-separated list (e.g., 39,128,256).",
        ),
        opts: str = typer.Option(
            None,
            "--opts",
            help="Override opts configuration as comma-separated list.",
        ),
        sector_opts: str = typer.Option(
            None,
            "--sector-opts",
            help="Override sector_opts configuration as comma-separated list.",
        ),
        foresight: str = typer.Option(
            None,
            "--foresight",
            help="Override foresight configuration (overnight, myopic, or perfect).",
        ),
        countries: str = typer.Option(
            None,
            "--countries",
            help="Override countries configuration as comma-separated list (e.g., DE,FR,ES).",
        ),
):
    """
    Run the PyPSA-Eur workflow with specified configuration.

    The workflow solves power system optimization problems based on the
    provided configuration parameters. Results are saved to the output directory.

    Examples:
    - pypsa-eur run --clusters 39 --foresight overnight
    - pypsa-eur run --configfile config/test/config.electricity.yaml --clusters 39,128
    - pypsa-eur run --name my_run --planning-horizons 2030,2040,2050
    """

    # Define a comprehensive mapping for rules to their human-readable actions
    rule_action_mapping = {
        # Data retrieval and preparation
        "retrieve_osm_prebuilt": "Retrieving OSM prebuilt network data",
        "retrieve_natura_raster": "Retrieving Natura 2000 raster data",
        "retrieve_cutout": "Retrieving weather cutout data",
        "retrieve_cost_data": "Retrieving technology cost data",
        "retrieve_ship_raster": "Retrieving shipping raster data",
        "retrieve_databundle": "Retrieving data bundle",
        "retrieve_egrigcs3_pm_production": "Retrieving EGRIGCS3 power production data",

        # Shape and geography building
        "build_shapes": "Building geographical shapes",
        "build_ship_raster": "Building shipping density raster",

        # Network building and preparation
        "base_network": "Building base network structure",
        "build_osm_network": "Building network from OpenStreetMap data",
        "clean_osm_network": "Cleaning OSM network topology",
        "build_bus_regions": "Building bus regions",
        "simplify_network": "Simplifying network topology",
        "cluster_network": "Clustering network nodes",
        "add_transmission_projects_and_dlr": "Adding transmission projects and dynamic line rating",
        "build_transmission_projects": "Building transmission projects",

        # Electricity system components
        "build_electricity_demand": "Building electricity demand time series",
        "build_electricity_demand_base": "Building base electricity demand",
        "build_powerplants": "Building power plant database",
        "build_line_rating": "Building dynamic line ratings",
        "add_electricity": "Adding electrical parameters",
        "prepare_network": "Preparing network for optimization",

        # Renewable energy
        "build_renewable_profiles": "Building renewable generation profiles",
        "determine_availability_matrix": "Determining renewable availability matrix",
        "build_renewable_potentials": "Building renewable energy potentials",
        "build_hydro_profile": "Building hydro generation profile",
        "build_temperature_profiles": "Building temperature profiles",
        "build_cop_profiles": "Building heat pump COP profiles",
        "build_solar_thermal_profiles": "Building solar thermal profiles",

        # Solving and optimization
        "solve_network": "Solving network optimization",
        "solve_operations_network": "Solving operational dispatch",
        "solve_elec_networks": "Solving electricity network",
        "solve_sector_network": "Solving sector-coupled network",

        # Sector coupling
        "prepare_sector_network": "Preparing sector-coupled network",
        "build_industrial_demand": "Building industrial demand",
        "build_heat_demand": "Building heat demand",
        "build_transport_demand": "Building transport demand",
        "build_biomass_potentials": "Building biomass potentials",
        "build_industry_sector": "Building industry sector model",
        "build_heating_sector": "Building heating sector model",
        "build_transport_sector": "Building transport sector model",

        # Post-processing and validation
        "make_summary": "Creating results summary",
        "plot_network": "Plotting network",
        "plot_summary": "Plotting summary statistics",
        "validate_network": "Validating network consistency",

        # Other rules
        "copy_config": "Copying configuration",
        "build_cutout": "Building weather cutout",
        "build_clustered_population_layouts": "Building clustered population layouts",
        "build_population_layouts": "Building population density layouts",
        "build_gas_network": "Building gas network",
        "build_gas_input_locations": "Building gas input locations",
        "build_industrial_production": "Building industrial production data",
        "build_industrial_energy_demand": "Building industrial energy demand",
        "build_district_heat_share": "Building district heating shares",
        "build_shipping_demand": "Building shipping demand",
        "build_aviation_demand": "Building aviation demand",
        "build_biomass_transport_costs": "Building biomass transport costs",
        "build_sequestration_potentials": "Building CO2 sequestration potentials",
        "build_co2_atmosphere": "Building atmospheric CO2 constraint",
        "build_salt_cavern_potentials": "Building salt cavern storage potentials",
    }

    # Use default config if not specified
    if not configfile:
        configfile = Path("config/config.default.yaml")
    else:
        # Ensure configfile exists
        if not configfile.exists():
            typer.secho(
                f"Error: Configuration file '{configfile}' does not exist.",
                fg=typer.colors.RED,
            )
            raise typer.Exit(code=1)

    # Generate job ID
    job_id = generate_job_id(configfile, name)
    typer.secho(f"Generated job ID: {job_id}", fg=typer.colors.CYAN)

    # Create a minimal config with only CLI overrides
    config_data = {}

    # Apply CLI overrides
    if name or job_id:
        config_data.setdefault("run", {})["name"] = job_id

    if snapshots:
        parts = snapshots.split(",")
        if len(parts) == 2:
            config_data.setdefault("snapshots", {})["start"] = parts[0].strip()
            config_data["snapshots"]["end"] = parts[1].strip()
        else:
            typer.secho(
                "Error: snapshots must be in format 'start,end'",
                fg=typer.colors.RED,
            )
            raise typer.Exit(code=1)

    if planning_horizons:
        ph_list = comma_separated_list(planning_horizons)
        config_data.setdefault("scenario", {})["planning_horizons"] = ph_list
        # Also set costs year to match the first planning horizon
        if ph_list:
            config_data.setdefault("costs", {})["year"] = ph_list[0]

    if clusters:
        config_data.setdefault("scenario", {})["clusters"] = comma_separated_list(clusters)

    if opts:
        config_data.setdefault("scenario", {})["opts"] = comma_separated_list(opts)

    if sector_opts:
        config_data.setdefault("scenario", {})["sector_opts"] = comma_separated_list(
            sector_opts
        )

    if foresight:
        if foresight not in ["overnight", "myopic", "perfect"]:
            typer.secho(
                f"Error: foresight must be 'overnight', 'myopic', or 'perfect', got '{foresight}'",
                fg=typer.colors.RED,
            )
            raise typer.Exit(code=1)
        config_data["foresight"] = foresight

    if countries:
        config_data["countries"] = comma_separated_list(countries)

    # Skip validation for now since we might be using EU-Flex schema
    # TODO: Add proper PyPSA-EUR schema validation
    typer.secho(
        "Skipping configuration validation (PyPSA-Eur schema not implemented)",
        fg=typer.colors.YELLOW,
    )

    typer.secho(
        "Initializing PyPSA-Eur workflow (interrupt anytime by pressing CTRL+C)",
        fg=typer.colors.GREEN,
    )

    # Write CLI-modified config to a temporary file
    cli_config_path = Path(tempfile.gettempdir()) / "config.pypsa_eur_cli.yaml"
    with open(cli_config_path, "w") as f:
        yaml.safe_dump(config_data, f)

    # Prepare configfile chain - CLI overrides should come LAST
    configfile_chain = ["--configfile"]

    # Start with user config or default config
    if configfile:
        configfile_chain.append(str(configfile.resolve()))
    else:
        default_config = Path("config/config.default.yaml").resolve()
        configfile_chain.append(str(default_config))

    # Add CLI overrides LAST so they take highest precedence
    configfile_chain.append(str(cli_config_path.resolve()))

    # Merge all configs to get final configuration
    final_config = {}
    for config_path in configfile_chain[1:]:  # Skip '--configfile'
        with open(config_path) as f:
            config_dict = yaml.safe_load(f)
            if config_dict:
                final_config = recursive_merge(final_config, config_dict)

    # Initialize the progress tracker with comprehensive rule mapping
    progress_tracker = ProgressTracker(rule_action_mapping=rule_action_mapping)

    try:
        project_root = Path.cwd()

        # Determine the target rule based on configuration
        # For electricity-only runs
        target_rule = "solve_elec_networks"

        # Build Snakemake command
        cmd = [
            "snakemake",
            target_rule,
            *configfile_chain,
            "-c",
            str(cores) if cores else "all",
            "--keep-incomplete",
        ]

        if verbose:
            typer.secho("Executing Snakemake command:", fg=typer.colors.BRIGHT_BLUE)
            typer.secho(f"  {' '.join(cmd)}", fg=typer.colors.BRIGHT_GREEN)
            typer.secho("\nWorkflow progress:", fg=typer.colors.CYAN)

            def run_verbose_snakemake():
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=str(project_root),
                )
                output_lines = []

                if process.stdout is None:
                    typer.secho(
                        "Error: Unable to get process output", fg=typer.colors.RED
                    )
                    return process.wait(), ""

                for line in process.stdout:
                    # Always print in verbose mode
                    tqdm.write(line.rstrip())
                    # Also parse for progress tracking
                    progress_tracker.parse_line(line)
                    output_lines.append(line)

                return_code = process.wait()
                progress_tracker.clean_up()
                return return_code, "".join(output_lines)

            # First attempt
            return_code, output = run_verbose_snakemake()

            # Check for lock and retry
            if "cannot be locked" in output or "LockException" in output:
                typer.secho(
                    "Detected Snakemake lock error. Attempting to unlock...",
                    fg=typer.colors.YELLOW,
                )
                unlock_cmd = [
                    "snakemake",
                    "--unlock",
                    "-c",
                    str(cores) if cores else "all",
                ]
                unlock_code = subprocess.call(unlock_cmd, cwd=str(project_root))

                if unlock_code == 0:
                    typer.secho(
                        "Successfully unlocked. Re-running Snakemake...",
                        fg=typer.colors.GREEN,
                    )
                    return_code, _ = run_verbose_snakemake()
                else:
                    typer.secho(
                        "Failed to unlock Snakemake directory.", fg=typer.colors.RED
                    )
                    raise typer.Exit(code=unlock_code)

            if return_code != 0:
                typer.secho("Workflow failed", fg=typer.colors.RED)
                raise typer.Exit(code=return_code)

            typer.secho("Workflow completed successfully", fg=typer.colors.GREEN)

        else:
            # Non-verbose mode - show progress with rule tracking
            typer.secho("Workflow progress:", fg=typer.colors.CYAN)

            process = subprocess.Popen(
                cmd,
                cwd=str(project_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )

            output_lines = []
            if process.stdout is not None:
                for line in process.stdout:
                    output_lines.append(line)
                    # Parse line for progress and rule tracking
                    progress_tracker.parse_line(line)

                process.wait()
                progress_tracker.clean_up()
                return_code = process.returncode

            output = "".join(output_lines)

            if "cannot be locked" in output or "LockException" in output:
                typer.secho(
                    "Detected Snakemake lock error. Attempting to unlock...",
                    fg=typer.colors.YELLOW,
                )
                unlock_cmd = ["snakemake", "--unlock"]
                unlock_code = subprocess.call(unlock_cmd, cwd=str(project_root))
                if unlock_code == 0:
                    typer.secho(
                        "Successfully unlocked. Re-running workflow...",
                        fg=typer.colors.GREEN,
                    )
                    return_code = subprocess.call(cmd, cwd=str(project_root))
                else:
                    typer.secho(
                        "Failed to unlock Snakemake directory.", fg=typer.colors.RED
                    )
                    raise typer.Exit(code=unlock_code)

            if return_code != 0:
                typer.secho("Workflow failed", fg=typer.colors.RED)
                raise typer.Exit(code=return_code)

            typer.secho("Workflow completed successfully", fg=typer.colors.GREEN)

    except FileNotFoundError as e:
        typer.secho(f"Error: {e}", fg=typer.colors.RED)
        typer.secho(
            "Hint: Is Snakemake installed and available?",
            fg=typer.colors.YELLOW,
        )
        raise typer.Exit(code=1)

    except KeyboardInterrupt:
        typer.secho("\nWorkflow interrupted by user", fg=typer.colors.YELLOW)
        raise typer.Exit(code=1)

    except Exception as e:
        typer.secho(f"Unexpected error: {e}", fg=typer.colors.RED)
        raise typer.Exit(code=1)


@app.callback()
def callback():
    """
    Main CLI entrypoint. Verifies that commands are run from the project root directory.
    """
    check_in_project_root()


if __name__ == "__main__":
    app()