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
    
    # Define a mapping for rules to their human-readable actions
    rule_action_mapping = {
        # Build electricity rules
        "build_shapes": "Building geographical shapes",
        "base_network": "Building base network",
        "build_bus_regions": "Building bus regions",
        "build_osm_network": "Building OSM network",
        "clean_osm_network": "Cleaning OSM network",
        "cluster_network": "Clustering network for",
        "simplify_network": "Simplifying network for",
        "add_electricity": "Adding electricity data for",
        "prepare_network": "Preparing network for",
        
        # Solve electricity rules
        "solve_network": "Solving network for",
        "solve_operations_network": "Solving operations network for",
        
        # Sector rules (if applicable)
        "prepare_sector_network": "Preparing sector network for",
        "solve_sector_network": "Solving sector network for",
        
        # Data retrieval
        "retrieve_osm_prebuilt": "Retrieving OSM prebuilt data",
        "retrieve_natura_raster": "Retrieving Natura raster data",
        "retrieve_cutout": "Retrieving cutout data",
        "retrieve_cost_data": "Retrieving cost data",
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
    
    # Load base configuration
    with configfile.open() as f:
        config_data = yaml.safe_load(f)
    
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
        config_data.setdefault("scenario", {})["planning_horizons"] = comma_separated_list(
            planning_horizons
        )
    
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
    
    # Validate configuration
    schema_path = "config/schema.yaml"
    
    if Path(schema_path).exists():
        try:
            typer.secho(
                "Validating PyPSA-Eur configuration", fg=typer.colors.YELLOW
            )
            validate(config_data, schema_path)
        except Exception as e:
            typer.secho(
                f"Configuration validation error: {e}",
                fg=typer.colors.RED,
            )
            # Continue anyway as schema might not be complete
            typer.secho(
                "Warning: Continuing despite validation error",
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
    
    # Prepare configfile chain
    default_config = Path("config/config.default.yaml").resolve()
    configfile_chain = ["--configfile"]
    
    # Always start with default config
    configfile_chain.append(str(default_config))
    
    # Add user config if different from default
    if configfile.resolve() != default_config:
        configfile_chain.append(str(configfile.resolve()))
    
    # Add CLI overrides
    configfile_chain.append(str(cli_config_path.resolve()))
    
    # Merge all configs to get final configuration
    final_config = {}
    for config_path in configfile_chain[1:]:  # Skip '--configfile'
        with open(config_path) as f:
            config_dict = yaml.safe_load(f)
            if config_dict:
                final_config = recursive_merge(final_config, config_dict)
    
    # Initialize the progress tracker
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
                    tqdm.write(line.rstrip())
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
            # Non-verbose mode
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