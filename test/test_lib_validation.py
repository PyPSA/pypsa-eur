from pathlib import Path

import pytest
import yaml
from jsonschema import ValidationError

from scripts.lib.validation import DeprecationError, validate_config


def test_validate_config_valid_configuration(tmp_path):
    """Test validation with a valid configuration."""
    # Create a temporary schema file
    schema = {
        "type": "object",
        "properties": {"name": {"type": "string"}, "age": {"type": "integer"}},
        "required": ["name", "age"],
    }

    schema_path = tmp_path / "schema.yaml"
    with open(schema_path, "w") as f:
        yaml.dump(schema, f)

    config = {"name": "John", "age": 30}

    validate_config(config, schema_path)


def test_validate_config_invalid_configuration(tmp_path):
    """Test validation with an invalid configuration."""
    schema = {
        "type": "object",
        "properties": {"name": {"type": "string"}, "age": {"type": "integer"}},
        "required": ["name", "age"],
    }

    schema_path = tmp_path / "schema.yaml"
    with open(schema_path, "w") as f:
        yaml.dump(schema, f)

    config = {"name": "John"}

    with pytest.raises(ValidationError):
        validate_config(config, schema_path)


def test_validate_config_wrong_type(tmp_path):
    """Test validation with wrong type in configuration."""
    schema = {
        "type": "object",
        "properties": {"name": {"type": "string"}, "age": {"type": "integer"}},
    }

    schema_path = tmp_path / "schema.yaml"
    with open(schema_path, "w") as f:
        yaml.dump(schema, f)

    config = {"name": "John", "age": "thirty"}

    with pytest.raises(ValidationError):
        validate_config(config, schema_path)


def test_validate_config_deprecated_keys(tmp_path):
    """Test validation with deprecated keys."""
    schema = {
        "type": "object",
        "properties": {
            "name": {"type": "string"},
            "age": {"type": "integer"},
            "old_field": {"type": "string", "deprecated": True},
        },
    }
    schema_path = tmp_path / "schema.yaml"
    with open(schema_path, "w") as f:
        yaml.dump(schema, f)

    config = {"name": "John", "age": 30, "old_field": "deprecated"}
    with pytest.raises(DeprecationError):
        validate_config(config, schema_path)


def test_validate_config_nonexistent_schema():
    """Test validation with non-existent schema file."""
    config = {"name": "John", "age": 30}

    with pytest.raises(FileNotFoundError):
        validate_config(config, "nonexistent_schema.yaml")


def test_validate_config_invalid_schema_format(tmp_path):
    """Test validation with invalid schema format."""
    schema_path = tmp_path / "invalid_schema.yaml"
    with open(schema_path, "w") as f:
        f.write("invalid: yaml: content")

    config = {"name": "John", "age": 30}

    with pytest.raises(yaml.YAMLError):
        validate_config(config, schema_path)


@pytest.mark.parametrize("schema_path", ["schema.yaml", Path("schema.yaml")])
def test_validate_config_path_types(tmp_path, schema_path):
    """Test validation with different path types."""
    schema = {"type": "object", "properties": {"name": {"type": "string"}}}

    full_path = tmp_path / schema_path
    with open(full_path, "w") as f:
        yaml.dump(schema, f)

    config = {"name": "John"}

    validate_config(config, full_path)
