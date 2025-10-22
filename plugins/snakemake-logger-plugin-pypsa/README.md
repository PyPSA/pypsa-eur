# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

# snakemake-logger-pypsa

A PyPSA logger plugin for Snakemake that provides detailed, human-readable output.

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
snakemake --logger pypsa
```

## Features

- Detailed, PyPSA-Eur tailored logging output

## Development

This package follows the Snakemake logger plugin interface defined in
`snakemake-interface-logger-plugins`.
