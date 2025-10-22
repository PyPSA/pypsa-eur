# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

# snakemake-logger-descriptive

A descriptive logger plugin for Snakemake that provides detailed, human-readable output.

## Installation

Install in development mode:

```bash
pip install -e .
```

Or with uv:

```bash
uv pip install -e .
```

## Usage

Use this logger plugin with Snakemake by specifying it via the `--logger` option:

```bash
snakemake --logger descriptive
```

## Features

- Detailed, descriptive logging output

## Development

This package follows the Snakemake logger plugin interface defined in
`snakemake-interface-logger-plugins`.
