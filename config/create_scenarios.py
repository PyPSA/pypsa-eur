# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import itertools

# Insert your config values that should be altered in the template.
template = """
scenario{scenario_number}:
    sector:
        carbon_: {config_value}

    config_section2:
        config_key2: {config_value2}
"""

# Define all possible combinations of config values.
# This must define all config values that are used in the template.
config_values = dict(config_values=["true", "false"], config_values2=[1, 2, 3, 4, 5])

combinations = [
    dict(zip(config_values.keys(), values))
    for values in itertools.product(*config_values.values())
]

# write the scenarios to a file
filename = "scenarios.yaml"
with open(filename, "w") as f:
    for i, config in enumerate(combinations):
        f.write(template.format(scenario_number=i, **config))
