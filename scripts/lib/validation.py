# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Validation function for input configuration files via JSON schema."""

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import yaml
from jsonschema import ValidationError
from jsonschema.validators import Draft202012Validator, extend


class DeprecationError(ValidationError):
    """Custom DeprecationError used during JSON schema validation."""

    pass


def _deprecated(
    validator: str, value: bool, instance: Any, schema: dict[str, Any]
) -> Iterator[DeprecationError]:
    yield DeprecationError(
        f"The key for value '{instance}' is deprecated: {schema['description']}",
        validator=validator,
        validator_value=value,
        instance=instance,
        schema=schema,
    )


CustomValidator = extend(Draft202012Validator, validators={"deprecated": _deprecated})


def validate_config(config: dict, schema_path: str | Path) -> None:
    """
    Validate the given configuration against the provided schema.

    Parameters
    ----------
    config : dict
        Configuration dictionary to validate.
    schema_path : str
        Path to the schema file.

    Raises
    ------
    ValidationError
        If the configuration is invalid according to the schema.
    DeprecationError
        If the configuration contains deprecated keys.

    """
    with open(schema_path) as schema_file:
        schema = yaml.safe_load(schema_file)

    custom_validator = CustomValidator(schema)
    custom_validator.validate(config)
