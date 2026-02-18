# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Config validation for PyPSA-EUR.

The schema is exported to both `config/config.default.yaml` and `config/schema.default.json`.
The json schema is also contributed to the schemastore.org and matches
`**/pypsa-eur*/config/*.yaml` to get IDE support without additional configuration.
"""

import re

from pydantic import ValidationError
from ruamel.yaml import YAML

from scripts.lib.validation.config._base import _registry
from scripts.lib.validation.config._schema import ConfigSchema


def validate_config(config: dict) -> ConfigSchema:
    """Validate config dict against schema."""
    config_schema = ConfigSchema
    name = config_schema._name.default
    for item in _registry:
        updater_config = item(config_schema)
        config_schema = updater_config.update()
        if updater_config.name:
            name += f".{updater_config.name}"
    validated_config = config_schema(**config)
    validated_config._name = name
    return validated_config


def generate_config_defaults(path: str = "config/config.{configname}.yaml") -> dict:
    """Generate config defaults YAML file and return the defaults dict."""
    from ruamel.yaml.comments import CommentedMap

    def convert_to_field_name(key: str) -> str:
        """Convert dash-case to snake_case for field lookup."""
        return key.replace("-", "_")

    # by_alias is needed to export dash-case instead of snake_case (which are some set aliases)
    # the goal should be to use snake_case consistently
    config = validate_config({})
    defaults = config.model_dump(by_alias=True)

    # Create YAML instance with custom settings
    yaml_writer = YAML()
    yaml_writer.version = (1, 1)  # Make sure to quote boolean-looking strings
    yaml_writer.default_flow_style = False
    yaml_writer.width = 4096  # Avoid line wrapping
    yaml_writer.indent(mapping=2, sequence=2, offset=0)

    # Custom string representer for controlling quote style
    def str_representer(dumper, data):
        """Use block style for multiline, quotes for special chars, plain otherwise."""
        TAG = "tag:yaml.org,2002:str"
        if "\n" in data:
            return dumper.represent_scalar(TAG, data, style="|")
        if data == "" or any(c in data for c in ":{}[]&*#?|-<>=!%@"):
            return dumper.represent_scalar(TAG, data, style='"')
        return dumper.represent_scalar(TAG, data, style="")

    yaml_writer.representer.add_representer(str, str_representer)

    # Create a CommentedMap to add comments
    data = CommentedMap()

    # Add yaml-language-server comment at the very top (before first key)
    data.yaml_set_start_comment(
        f"yaml-language-server: $schema=./schema.{config._name}.json"
    )

    for key, value in defaults.items():
        data[key] = value

        field_name = convert_to_field_name(key)
        docs_url = f"https://pypsa-eur.readthedocs.io/en/latest/configuration.html#{field_name}"
        data.yaml_set_comment_before_after_key(key, before=f"\ndocs in {docs_url}")

    # Write to file
    with open(path.format(configname=config._name), "w") as f:
        yaml_writer.dump(data, f)

    return defaults


def generate_config_schema(path: str = "config/schema.{configname}.json") -> dict:
    """Generate JSON schema file and return the schema dict."""
    import json
    import math

    def resolve_refs(obj: dict, defs: dict) -> dict:
        """Resolve nested schema references to show them nicely in the documentation."""
        if isinstance(obj, dict):
            if "$ref" in obj:
                ref_path = obj["$ref"]  # "#/$defs/RunConfig
                ref_name = ref_path.split("/")[-1]
                if ref_name in defs:
                    resolved = resolve_refs(defs[ref_name].copy(), defs)
                    # Keep description from the reference
                    if "description" in obj and "description" not in resolved:
                        resolved["description"] = obj["description"]
                    return resolved
            return {k: resolve_refs(v, defs) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [resolve_refs(item, defs) for item in obj]
        return obj

    def sanitize_for_json(obj):
        """Replace infinity values with None for valid JSON."""
        if isinstance(obj, dict):
            return {k: sanitize_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [sanitize_for_json(v) for v in obj]
        elif isinstance(obj, float) and math.isinf(obj):
            return None
        return obj

    def remove_nested_titles(obj, is_root=True):
        """Remove nested titles (e.g. model class names)."""
        if isinstance(obj, dict):
            result = {}
            for k, v in obj.items():
                if k == "title" and not is_root:
                    continue
                result[k] = remove_nested_titles(v, is_root=False)
            return result
        elif isinstance(obj, list):
            return [remove_nested_titles(item, is_root=False) for item in obj]
        return obj

    def remove_object_type(obj, is_root=True):
        """Remove 'type: object' from nested objects (redundant when properties exist)."""
        if isinstance(obj, dict):
            result = {}
            for k, v in obj.items():
                if (
                    k == "type"
                    and v == "object"
                    and not is_root
                    and "properties" in obj
                ):
                    continue
                result[k] = remove_object_type(v, is_root=False)
            return result
        elif isinstance(obj, list):
            return [remove_object_type(item, is_root=False) for item in obj]
        return obj

    def convert_rst_to_markdown(obj):
        """Convert RST-style links in 'description' to Markdown in 'markdownDescription'."""

        def rst_to_md(text):
            """Convert RST link format `Link Text <URL>`_ to Markdown [Link Text](URL)."""
            # Pattern matches: `Link Text <URL>`_
            pattern = r"`([^<>`]+)\s*<([^>]+)>`_"
            return re.sub(pattern, r"[\1](\2)", text)

        if isinstance(obj, dict):
            result = {}
            for k, v in obj.items():
                if k == "description" and isinstance(v, str) and "`" in v and "<" in v:
                    result[k] = v
                    md_text = rst_to_md(v)
                    if md_text != v:
                        result["markdownDescription"] = md_text
                else:
                    result[k] = convert_rst_to_markdown(v)
            return result
        elif isinstance(obj, list):
            return [convert_rst_to_markdown(item) for item in obj]
        return obj

    config = validate_config({})
    schema = config.model_json_schema()
    defs = schema.get("$defs", {})
    schema = resolve_refs(schema, defs)
    schema = sanitize_for_json(schema)
    schema = remove_nested_titles(schema)
    schema = remove_object_type(schema)
    schema = convert_rst_to_markdown(schema)
    with open(path.format(configname=config._name), "w") as f:
        json.dump(schema, f, indent=2)
        f.write("\n")
    return schema


__all__ = [
    "ConfigSchema",
    "validate_config",
    "generate_config_defaults",
    "generate_config_schema",
    "ValidationError",
]
