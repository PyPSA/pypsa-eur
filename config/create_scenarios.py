# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# This script helps to generate a scenarios.yaml file for PyPSA-Eur.
# You can modify the template to your needs and define all possible combinations of config values that should be considered.

if "snakemake" in globals():
    filename = snakemake.output[0]  # noqa: F821
else:
    filename = "../config/scenarios.yaml"

import itertools

# Insert your config values that should be altered in the template.
# Change `config_section` and `config_section2` to the actual config sections.
template = """
scenario{scenario_number}:
    config_section:
        config_key: {config_value}
    config_section2:
        config_key2: {config_value2}
"""

# Define all possible combinations of config values.
# This must define all config values that are used in the template.
config_values = dict(config_value=["true", "false"], config_value2=[1, 2, 3, 4])

combinations = [
    dict(zip(config_values.keys(), values))
    for values in itertools.product(*config_values.values())
]

with open(filename, "w") as f:
    for i, config in enumerate(combinations):
        f.write(template.format(scenario_number=i, **config))
