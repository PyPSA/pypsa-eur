# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
This script extracts the results of the GSA from the model runs and calculates the GSA metrics.

Please see the file config/GSA.yaml for an explanation of the workflow.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from SALib.analyze import morris as analyze_morris
from SALib.plotting import morris as plot_morris


def create_salib_problem(parameters: list) -> dict:
    """
    Load the GSA configuration from `config/GSA.yaml`.

    Returns
    -------
    dict
        GSA configuration dictionary.
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


def extract_results():
    """
    Extract results from model runs and save them as CSV files for GSA analysis.
    """
    resultsfolder = Path("results")
    gsa_config = get_gsa_config()
    result_variables = gsa_config.get("results", [])
    results_dir = Path("GSA/results")
    results_dir.mkdir(parents=True, exist_ok=True)
    for variable, params in result_variables.items():
        variable_dict = {}
        for folder in resultsfolder.iterdir():
            if "modelrun" in folder.name:
                run = folder.name
                for folder in folder.iterdir():
                    if folder.is_dir() and folder.name == "csvs":
                        for file in folder.iterdir():
                            if file.is_file() and file.name == params.get("csv").get(
                                "filename"
                            ):  # filters the csv file of the results variable
                                identification_string = params.get("csv").get(
                                    "identification_column_entries"
                                )
                                identification_entries = identification_string.split(
                                    ","
                                )
                                csv = file.open()  # opens the csv file
                                # all the len(identification_entries) columns need to match the identification_entries
                                # the last column is the one that is being extracted
                                for line in csv:
                                    columns = line.split(",")
                                    if all(
                                        [
                                            columns[i] == identification_entries[i]
                                            for i in range(len(identification_entries))
                                        ]
                                    ):
                                        if columns[-1] == "nan":
                                            value = 0
                                        else:
                                            value = float(columns[-1])
                                        multiplier = float(params.get("multiplier"))
                                        variable_dict[run] = value * multiplier
                                        break
                                csv.close()

        variable_dict = dict(
            sorted(
                variable_dict.items(),
                key=lambda item: int("".join(filter(str.isdigit, item[0]))),
            )
        )  # sorts the dictionary by the run number
        output_file = results_dir / f"{variable}.csv"
        with open(output_file, "w") as f:
            f.write("run,value\n")
            for run, value in variable_dict.items():
                f.write(f"{run},{value}\n")
        f.close()


def calculate_GSA_metrics():
    """
    Calculate GSA metrics and generate sensitivity plots.
    """
    gsa_config = get_gsa_config()
    parameters = gsa_config["parameters"]
    sample = "GSA/morris_sample.txt"

    # create directory for the GSA results
    Path("GSA/SA_results").mkdir(parents=True, exist_ok=True)
    Path("GSA/SA_plots").mkdir(parents=True, exist_ok=True)

    X = np.loadtxt(sample, delimiter=",")
    problem = create_salib_problem(parameters)
    for variable, params in gsa_config.get("results", {}).items():
        results = pd.read_csv(f"GSA/results/{variable}.csv")
        Y = results["value"].values
        if len(X) != len(Y):
            raise ValueError(
                f"Length of sample ({len(X)}) and results ({len(Y)}) do not match. Try changing the identification_column_entries in the GSA.yaml file."
            )
        print(f"GSA for {variable}: ")
        Si = analyze_morris.analyze(problem, X, Y, print_to_console=True)
        print("\n")
        Si.to_df().to_csv(f"GSA/SA_results/{variable}_SA_results.csv")

        title = variable.capitalize()
        unit = params.get("unit", "")

        # Check if covariance plot should be included
        if gsa_config["general"].get("covariance_plot", True):
            fig, axs = plt.subplots(2, figsize=(20, 20))
            plot_morris.horizontal_bar_plot(axs[0], Si, unit=unit)
            plot_morris.covariance_plot(axs[1], Si, unit=unit)
            axs[1].set_xlabel(f"µ in {unit}", fontsize=18)
            axs[1].tick_params(axis="both", which="major", labelsize=14)
        else:
            fig, axs = plt.subplots(1, figsize=(16, 10))
            plot_morris.horizontal_bar_plot(axs, Si, unit=unit)
            axs.set_xlabel(f"µ* in {unit}", fontsize=18)
            axs.tick_params(axis="both", which="major", labelsize=14)
        fig.suptitle(title, fontsize=24)

        fig.savefig(f"GSA/SA_plots/{variable}.png", bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    extract_results()
    calculate_GSA_metrics()
