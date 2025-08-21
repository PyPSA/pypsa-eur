import hashlib
import json
import time
from datetime import datetime
from pathlib import Path

import typer
import yaml


def generate_job_id(configfile: Path | None = None, name: str | None = None) -> str:
    """Generate a unique job ID based on the config and the current timestamp."""
    if name:
        # If name is provided via CLI, use it directly
        return name
    
    if configfile and configfile.exists():
        # Load the configuration file
        with open(configfile) as f:
            config_data = yaml.safe_load(f)
        
        # Extract the original run name from the config
        original_run_name = config_data.get("run", {}).get("name", "pypsa_eur")
    else:
        original_run_name = "pypsa_eur"
    
    # Generate timestamp-based unique ID
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create a short hash for uniqueness
    hash_input = f"{original_run_name}{timestamp}{time.time()}"
    short_hash = hashlib.sha256(hash_input.encode()).hexdigest()[:6]
    
    # Combine everything
    final_job_id = f"{timestamp}_{original_run_name}_{short_hash}"
    
    return final_job_id


def comma_separated_list(value: str) -> list:
    """Parse comma-separated values into a list."""
    if isinstance(value, str):
        # Try to parse as integers first
        try:
            return [int(v.strip()) for v in value.split(",") if v.strip()]
        except ValueError:
            # If not all integers, return as strings
            return [v.strip() for v in value.split(",") if v.strip()]
    elif isinstance(value, (list, tuple)):
        return list(value)
    elif isinstance(value, (int, float)):
        return [value]
    else:
        raise ValueError(
            f"Invalid input: expected comma-separated values, got '{value}'"
        )


def check_in_project_root():
    """Raise an error if not in project root directory."""
    expected_file = Path("config/config.default.yaml")
    if not expected_file.exists():
        typer.secho(
            f"Error: Expected to find {expected_file} in current directory. "
            "Make sure you are running the CLI from the root of the PyPSA-Eur project.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)


def recursive_merge(d1: dict, d2: dict) -> dict:
    """Recursively merge dictionaries, values in d2 overwrite those in d1."""
    for key, value in d2.items():
        if isinstance(value, dict) and key in d1 and isinstance(d1[key], dict):
            recursive_merge(d1[key], value)  # Recursively merge dictionaries
        else:
            d1[key] = value  # Overwrite or add the value
    return d1