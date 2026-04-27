# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for config schema synchronization."""

import difflib
import tempfile
from pathlib import Path

import pytest

from scripts.lib.validation.config import (
    generate_config_defaults,
    generate_config_schema,
    validate_config,
)


@pytest.fixture(scope="module")
def config_file():
    config = validate_config({})
    config_filename = Path(f"config/config.{config._name}.yaml")
    return config_filename


@pytest.fixture(scope="module")
def schema_file():
    config = validate_config({})
    schema_filename = Path(f"config/schema.{config._name}.json")
    return schema_filename


def _check_file_in_sync(existing_path: Path, generate_func, file_type: str):
    """Check if a file is in sync with the generated content."""
    assert existing_path.exists(), f"{existing_path} does not exist"

    existing_content = existing_path.read_text()

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=f".{file_type}", delete=False
    ) as f:
        temp_path = f.name

    try:
        generate_func(temp_path)
        generated_content = Path(temp_path).read_text()
    finally:
        Path(temp_path).unlink(missing_ok=True)

    if existing_content != generated_content:
        diff = difflib.unified_diff(
            existing_content.splitlines(keepends=True),
            generated_content.splitlines(keepends=True),
            fromfile=f"{existing_path} (existing)",
            tofile=f"{existing_path} (generated)",
        )
        diff_str = "".join(diff)
        print(f"\n{'=' * 80}")
        print(f"DIFF: {existing_path}")
        print("=" * 80)
        print(diff_str)
        print("=" * 80)
        pytest.fail(
            f"{existing_path} is out of sync with the Pydantic schema. "
            "Run 'pixi run generate-config' to update. See diff above."
        )


def test_config_default_yaml_in_sync(config_file):
    """Test that config/config.default.yaml is in sync with Pydantic schema."""

    _check_file_in_sync(
        config_file,
        generate_config_defaults,
        "yaml",
    )


def test_config_schema_json_in_sync(schema_file):
    """Test that config/schema.default.json is in sync with Pydantic schema."""

    _check_file_in_sync(
        schema_file,
        generate_config_schema,
        "json",
    )
