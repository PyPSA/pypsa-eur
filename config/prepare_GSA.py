# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
This script is used to generate the scenarios for the Global Sensitivity Analysis (GSA) using the SALib library.
The script reads the parameters from the GSA.yaml file and generates the scenarios for the GSA.
The scenarios are saved in the config/GSA_runs.yaml file.

See the file config/GSA.yaml for an explanation of the workflow.
"""

from pathlib import Path

import numpy as np
import yaml
from SALib.sample import morris


def create_salib_problem(parameters: list) -> dict:
    """
    Create a SALib problem definition from the parameter configuration.

    Parameters
    ----------
    parameters : dict
        Dictionary of parameters with their ranges and group names.

    Returns
    -------
    dict
        SALib problem definition with names, bounds, and groups.
    """
    problem = {}
    problem["num_vars"] = len(
        parameters
    )  # this is the number of parameters (to be extracted of the configfile)
    if problem["num_vars"] <= 1:
        raise ValueError(
            f"Must define at least two variables in problem. User defined "
            f"{problem['num_vars']} variable(s)."
        )

    names = []
    bounds = []
    groups = []
    for param_name, param_details in parameters.items():
        names.append(param_name)  # Use the parameter name directly
        groups.append(param_details["groupname"])
        min_value = param_details["min"]
        max_value = param_details["max"]
        bounds.append([min_value, max_value])

    problem["names"] = names
    problem["bounds"] = bounds
    problem["groups"] = groups
    num_groups = len(set(groups))
    if num_groups <= 1:
        raise ValueError(
            f"Must define at least two groups in problem. User defined "
            f"{num_groups} group(s)."
        )
    return problem


def set_nested_value(d, keys, value):
    for key in keys[:-1]:
        d = d.setdefault(key, {})
    d[keys[-1]] = value


def extract_keys(d, keys):
    for k, v in d.items():
        keys.append(k)
        if isinstance(v, dict):
            extract_keys(v, keys)


def generate_scenarios(sample, output_file: str, parameters: dict):
    """
    Generate GSA scenarios from the sampled parameter values and save them to a YAML file.

    Parameters
    ----------
    sample : np.ndarray
        Array of sampled parameter values.
    output_file : str
        Path to the output YAML file.
    parameters : dict
        Dictionary of parameter definitions.
    """
    scenarios = {}
    for i, row in enumerate(sample):
        scenario = {}
        for j, (param_name, param_details) in enumerate(parameters.items()):
            value = float(row[j])
            config_file_location = param_details["config_file_location"]

            # Convert the nested dictionary structure to a flat list of keys
            keys = []
            extract_keys(config_file_location, keys)

            # Set the value directly in the scenario dictionary without the parameter name level
            set_nested_value(scenario, keys, value)

        scenarios[f"modelrun_{i}"] = scenario

    yaml_content = ""
    for i in range(len(scenarios)):
        key = f"modelrun_{i}"
        if key in scenarios:
            scenario_yaml = yaml.dump({key: scenarios[key]}, default_flow_style=False)
            yaml_content += scenario_yaml + "\n"

    with open(output_file, "w") as f:
        f.write(yaml_content)


def get_gsa_config() -> dict:
    """
    Load the GSA configuration from `config/GSA.yaml`.

    Returns
    -------
    dict
        GSA configuration dictionary.
    """
    config_path = Path("config/GSA.yaml")
    if not config_path.exists():
        config_path = Path("config/GSA.default.yaml")
        if not config_path.exists():
            raise FileNotFoundError(
                "No GSA configuration file found. Please create a GSA.yaml file in the config directory."
            )
        print(
            f"File config/GSA.yaml not found. Using default configuration ({config_path})."
        )
    with config_path.open() as f:
        return yaml.safe_load(f)


def main():
    gsa_config = get_gsa_config()
    sample_file = "GSA/morris_sample.txt"

    replicates = gsa_config["general"]["replicates"]
    parameters = gsa_config["parameters"]

    problem = create_salib_problem(parameters)
    sample = morris.sample(
        problem,
        N=100,
        optimal_trajectories=replicates,
        local_optimization=True,
        seed=42,
    )

    Path("GSA").mkdir(parents=True, exist_ok=True)
    np.savetxt(sample_file, sample, delimiter=",")

    parameters = gsa_config["parameters"]
    output_file = "config/GSA_runs.yaml"

    morris_sample = np.loadtxt(sample_file, delimiter=",")
    generate_scenarios(morris_sample, output_file, parameters)


if __name__ == "__main__":
    main()
